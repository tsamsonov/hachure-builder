/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.hachurebuilder;

import autolab.grid.ESRIASCIIGridHeader;
import autolab.grid.HachureBuilder;
import com.vividsolutions.jts.geom.LineString;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.Transaction;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.FeatureCollection;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import autolab.grid.EsriASCIIGridReader;
import com.vividsolutions.jts.awt.ShapeReader;
import com.vividsolutions.jts.geom.CoordinateSequence;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.MultiLineString;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.concurrent.ExecutionException;
import marchingsquares.Algorithm;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTSFactoryFinder;
/**
 *
 * @author tsamsonov
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        args = new String[11];
                
        args[0] = "/Volumes/Data/Work/__RGO/Hachures/kondyor.asc";
        args[1] = "/Volumes/Data/Work/__RGO/Hachures/kondyor.shp";
        args[2] = "0";
        args[3] = "100";
        args[4] = "10"; 
        args[5] = "5";
        args[6] = "1";
        args[7] = "2"; 
        args[8] = "40";
        args[9] = "10";
        args[10] = "0";
                
        int nargs = args.length;
        if(nargs != 11){
            System.err.println("ARGUMENTS ERROR: " + nargs + " arguments given. 11 expected");
        } else {
            try{
                String dempath = args[0];
                String outpath = args[1];
                double z0 = Double.parseDouble(args[2]);
                double dz = Double.parseDouble(args[3]);
                double dist = Double.parseDouble(args[4]);
                double minDist = Double.parseDouble(args[5]);
                double step = Double.parseDouble(args[6]);
                double minSlope = Double.parseDouble(args[7]);
                double maxTurn = Double.parseDouble(args[8]);
                int maxDepth = Integer.parseInt(args[9]);
                boolean isVariable = Boolean.parseBoolean(args[10]);
                
                // READ INPUT
                
                File file = new File( dempath );

//                AbstractGridFormat format = GridFormatFinder.findFormat( file );
//                GridCoverage2DReader reader = format.getReader( file );
                
//                GeoTiffReader reader = new GeoTiffReader(file);
//                GridCoverage2D dem = (GridCoverage2D) reader.read(null);
                
                double[][] mtx = EsriASCIIGridReader.readDoubleMatrix(file);
                
                FileInputStream fis = new FileInputStream(file.getAbsolutePath());
                Scanner scanner = EsriASCIIGridReader.createUSScanner(fis);       
                ESRIASCIIGridHeader header = new ESRIASCIIGridHeader();
                header.readHeader(scanner);
                
                double lon0 = header.west;
                double lat0 = header.south;
                int cellsX = header.cols;
                int cellsY = header.rows;
                double lon1 = lon0 + cellsX * header.cellSize;
                double lat1 = lat0 + cellsY * header.cellSize;
                
                double[] levels = new double[4]; // create 1-D array of thresholds.
                levels[0] = 600;
                levels[1] = 800;
                levels[2] = 1000;
                levels[3] = 1200;
                    
                Algorithm alg = new Algorithm();
                GeneralPath[] isolines = alg.buildContours(mtx, levels);
                
                AffineTransform xf = new AffineTransform();
       
                xf.translate(lon0, lat0);
                xf.scale((lon1 - lon0) / (cellsX - 1), (lat1 - lat0) / (cellsY - 1));
                xf.translate(-1, -1); // Because MxN data was padded to (M+2)x(N+2).
                for (int i = 0; i < isolines.length; i++) {
                    isolines[i].transform(xf); // Permanent mapping to world coords.
                }
//                SimpleFeatureTypeBuilder typeBuilder = new SimpleFeatureTypeBuilder();        
//                typeBuilder.setName("Contours");                
//                Class geomType = LineString.class;
//                typeBuilder.add("Location", geomType);
//                typeBuilder.add("Level", Double.class);
//                final SimpleFeatureType TYPE = typeBuilder.buildFeatureType();
                GeometryFactory gFactory = JTSFactoryFinder.getGeometryFactory();
//                SimpleFeatureBuilder fBuilder = new SimpleFeatureBuilder(TYPE);
                
                ArrayList<Geometry> contours = new ArrayList<>();
                for(int i = 0; i < isolines.length; i++){
                    Geometry cnt = ShapeReader.read(isolines[i].getPathIterator(null), gFactory);
                    contours.add(cnt);
                }

//                
                GridCoverage2D dem = EsriASCIIGridReader.read(file);
                
//                GridCoverageReader reader = new ArcGridReader(file);
//                GridCoverage2D dem = (GridCoverage2D) reader.read(null);
//                CoordinateReferenceSystem crs = dem.getCoordinateReferenceSystem2D();
                
                // CONSTRUCT HACHURES
                
                HachureBuilder builder = new HachureBuilder(dem, z0, dz);
                
                FeatureCollection hachures = builder.getHachures(dist, minDist, 
                        step, minSlope, maxTurn, maxDepth, isVariable);
                
                // WRITE OUTPUT
                
                File outfile = new File( outpath );
                ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();
                Map<String, Serializable> params = new HashMap<>();
                params.put("url", outfile.toURI().toURL());
                params.put("create spatial index", Boolean.TRUE);

                ShapefileDataStore dataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);
                                        
                final SimpleFeatureType TYPE = (SimpleFeatureType)hachures.getSchema();
                
//                SimpleFeatureTypeBuilder sftb = new SimpleFeatureTypeBuilder();
//                sftb.init(TYPE);
//                sftb.setCRS(crs);
//                
//                dataStore.createSchema(sftb.buildFeatureType());
               
                dataStore.createSchema(TYPE);
                String typeName = dataStore.getTypeNames()[0];
                SimpleFeatureSource featureSource = dataStore.getFeatureSource(typeName);
                SimpleFeatureStore featureStore = (SimpleFeatureStore) featureSource;

                SimpleFeatureIterator features = (SimpleFeatureIterator) hachures.features();
                while(features.hasNext()){
                    SimpleFeature f = features.next();
                    LineString l = (LineString)f.getDefaultGeometry();
                    if(l.isClosed()){
                        int c = 0;
                    }
                    if(l.getStartPoint().getCoordinate().x == l.getEndPoint().getCoordinate().x){
                        int c = 0;
                    }
                }

                Transaction transaction = new DefaultTransaction("create");
                try{
                    featureStore.setTransaction(transaction);
                    featureStore.addFeatures(hachures);
                    transaction.commit();
                } catch (IOException problem){
                    transaction.rollback();
                } finally {
                    transaction.close();
                }
                    
            } catch (MalformedURLException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            } catch(NumberFormatException e){
                System.err.println("ARGUMENTS ERROR: Invalid type");
            } catch(IOException e){
                System.err.println(e.getMessage());
                e.printStackTrace();
                System.err.println("READ ERROR: Failed to read DEM file");
            } catch (InterruptedException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            } catch (ExecutionException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
}
