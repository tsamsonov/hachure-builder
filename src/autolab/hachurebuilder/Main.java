/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.hachurebuilder;

import autolab.grid.HachureBuilder;
import java.io.File;
import java.io.IOException;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.io.AbstractGridFormat;
import org.geotools.coverage.grid.io.GridCoverage2DReader;
import org.geotools.coverage.grid.io.GridFormatFinder;
import org.geotools.feature.FeatureCollection;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

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
        int nargs = args.length;
        if(nargs != 10){
            System.err.println("ARGUMENTS ERROR: " + nargs + " arguments given. 10 expected");
        } else {
            try{
                String demfile = args[0];
                String output = args[1];
                double z0 = Double.parseDouble(args[2]);
                double dz = Double.parseDouble(args[3]);
                double dist = Double.parseDouble(args[4]);
                double minDist = Double.parseDouble(args[5]);
                double step = Double.parseDouble(args[6]);
                double minSlope = Double.parseDouble(args[7]);
                double maxTurn = Double.parseDouble(args[8]);
                int maxDepth = Integer.parseInt(args[9]);
                boolean isVariable = Boolean.parseBoolean(args[10]);
                
                File file = new File( demfile );

                AbstractGridFormat format = GridFormatFinder.findFormat( file );
                GridCoverage2DReader reader = format.getReader( file );
                
                GridCoverage2D dem = (GridCoverage2D) reader.read(null);
                CoordinateReferenceSystem crs = dem.getCoordinateReferenceSystem2D();
                
                HachureBuilder builder = new HachureBuilder(dem, z0, dz);
                
                FeatureCollection hachures = builder.getHachures(dist, minDist, 
                        step, minSlope, maxTurn, maxDepth, isVariable);
                
                
                
                
            } catch(NumberFormatException e){
                System.err.println("ARGUMENTS ERROR: Invalid type");
            } catch(IOException e){
                System.err.println("READ ERROR: Failed to read DEM file");
            }
        }
    }
    
}
