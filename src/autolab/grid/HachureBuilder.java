/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.grid;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateSequence;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.factory.CommonFactoryFinder;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.process.raster.ContourProcess;

import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.filter.Filter;
import org.opengis.filter.FilterFactory2;

import autolab.math.AutolabMath;
import java.lang.Math;
import java.util.Arrays;

/**
 *
 * @author tsamsonov
 */
public class HachureBuilder {
    
    double Z[][];
    double minSlope;
    double maxTurn;
    double dist; // hachure distances
    double mindist; // hachure min distance (for vectors)
    double step; // approximation step
    int maxDepth;
    String HField; // the field containing height of the contour
    boolean isVector; // streams are extracted from vector field grids
    
    GridProcessor gp; // scalar field
    GridProcessor gpx; // dx component of vector field
    GridProcessor gpy; // dy component of vector field
    
    ArrayList<CoordinateSequence> hachuresDown;
    ArrayList<CoordinateSequence> hachuresUp;
    ArrayList<Coordinate> seeds;
    SimpleFeatureCollection contours;
    
    // for variable step hachures
    ArrayList<Double> distances;
    
    // contour levels;
    ArrayList<Double> levels;
    
    
    ArrayList<Double> hachureSpeeds;
    
    ArrayList<Double> hachureLevels;
    
    /**
     * Hachures constructor
     * @param dem elevation data
     * @param z0 start elevation level
     * @param dz contour interval
     */
    public HachureBuilder(GridCoverage2D dem,
                          Double z0,
                          Double dz){
        this.gp = new GridProcessor(dem);
        this.isVector = false;
        
        Double zmin = gp.getHeader().zmin;
        Double zmax = gp.getHeader().zmax;
        
        // Define the basic level
        double level = zmin + (z0 - zmin) % dz + dz * ((zmin > z0) ? 1 : 0);
        
        int nlevels = (int)Math.floor((zmax-level)/z0) + 1;
        double[] contourLevels = new double[nlevels];
        levels =  new ArrayList<>(nlevels);
        
        for(int i = 0; i < nlevels; i++){
            contourLevels[i] = level;
            levels.set(i, level);
            level += dz;
        }
        
        this.contours = ContourProcess.process(dem, 0, contourLevels, dz, true, true, null, null);
        this.HField = "value";
        
        // Sort levels in decreasing order
        Comparator cmp = Collections.reverseOrder();
        Collections.sort(levels,cmp);
        
        // add minimum and maximum
        levels.add(zmin);
        levels.add(0,zmax);
        
    }
    
    /**
     * Hachures constructor
     * @param dem Digital elevation model
     * @param contours contours (sorted high to low)
     * @param HField height field
     */
    public HachureBuilder(GridCoverage2D dem, 
                    FeatureCollection contours,
                    String HField)
    {
        this.gp = new GridProcessor(dem);
        this.contours = (SimpleFeatureCollection)contours;
        this.HField = HField;
        this.isVector = false;
        
        // Get contour levels values
        
        levels = new ArrayList<>();
        SimpleFeatureIterator contIterator = (SimpleFeatureIterator)contours.features();
        while(contIterator.hasNext()){
            SimpleFeature cont = contIterator.next();
            Double level = (Double)cont.getAttribute(HField);
            
            if(!levels.contains(level))
                levels.add(level);
        }
        
        // Sort levels in decreasing order
        Comparator cmp = Collections.reverseOrder();
        Collections.sort(levels,cmp);
        
        // add minimum and maximum values
        levels.add(0, gp.getHeader().zmax);
        levels.add(gp.getHeader().zmin);
    }
    
    /**
     * Hachures constructor
     * @param dem Digital elevation model
     * @param dxdem
     * @param dydem
     * @param isVariable
     * @param contours contours (sorted high to low)
     * @param HField height field
     */
    public HachureBuilder(GridCoverage2D dem,
                          GridCoverage2D dxdem,
                          GridCoverage2D dydem,
                    FeatureCollection contours,
                    String HField)
    {
        gp = new GridProcessor(dem);
        gpx = new GridProcessor(dxdem);
        gpy = new GridProcessor(dydem);
        this.contours = (SimpleFeatureCollection)contours;
        this.HField = HField;
        isVector = true;
        
        // Get contour levels values
        
        levels = new ArrayList<>();
        SimpleFeatureIterator contIterator = (SimpleFeatureIterator)contours.features();
        while(contIterator.hasNext()){
            SimpleFeature cont = contIterator.next();
            Double level = (Double)cont.getAttribute(HField);
            
            if(!levels.contains(level))
                levels.add(level);
        }
        
        // Sort levels in decreasing order
        Comparator cmp = Collections.reverseOrder();
        Collections.sort(levels,cmp);
        
        // add minimum and maximum values
        levels.add(0, gp.getHeader().zmax);
        levels.add(gp.getHeader().zmin);
        
    }

    /**
     * Calculates hachures and writes them into FeatureCollection
     * @param dist
     * @param step
     * @param minSlope
     * @param maxTurn
     * @param maxDepth
     * @return 
     */
    public FeatureCollection getHachures(double dist, 
                                         double minDist,
                                         double step,
                                         double minSlope,
                                         double maxTurn,
                                         int maxDepth,
                                         boolean isVariable){
        hachuresDown = new ArrayList<>();
        hachuresUp = new ArrayList<>();
        hachureLevels = new ArrayList<>();
        
        seeds = new ArrayList<>();
        
        this.minSlope = minSlope*Math.PI/180;
        this.maxTurn = maxTurn*Math.PI/180;
        this.dist = dist*gp.getHeader().res; // hachure distance
        this.mindist = minDist*gp.getHeader().res; // hachure min distance
        this.step = step*gp.getHeader().res; // approximation step
        this.maxDepth = maxDepth;

        // prepare distances
        distances = new ArrayList<>();
        double delta = isVariable ? (dist-minDist)/levels.size() : 0;
        for(int i = 0; i < levels.size()-2; i++){
            distances.add((minDist + delta*i)*gp.getHeader().res);
        }
        
        // process each contour level
        FilterFactory2 ff = CommonFactoryFinder.getFilterFactory2();
        
        int nFirstPrev = 0; // first index of previously constructed (upper) level
        int n = 0; //number of currently constructed hachures
        for(int i = 1; i < levels.size() - 1; i++){
            double maxZ = levels.get(i-1);
            double curZ = levels.get(i);
            double minZ = levels.get(i+1);
            
            Filter filter = ff.equals(ff.property(HField), ff.literal(curZ));
            
            SimpleFeatureCollection currLevel = contours.subCollection(filter);
            
            int nf = currLevel.size();
            
            SimpleFeatureIterator levelIterator = 
                    (SimpleFeatureIterator)currLevel.features();
            
            int nFirst  = hachuresDown.size();
            
            double curDist = distances.get(i-1);
            
            System.out.print("Processing level " + i + " out of " + levels.size());
            while(levelIterator.hasNext()){
                SimpleFeature line = levelIterator.next();
                LineString cont = (LineString)(((MultiLineString)line.getDefaultGeometry()).getGeometryN(0));
                
                int nFirstCont = hachuresDown.size(); // first index of currently constructed level
                // get main hachures
                getHachures(cont, curDist, maxZ, curZ, minZ, nFirstPrev);
                
                // insert additional hachures
                insertHachures(curDist, minZ, curZ, maxZ, nFirstCont, nFirst);
            }
            nFirstPrev = nFirst;
        }
        
        SimpleFeatureTypeBuilder typeBuilder = new SimpleFeatureTypeBuilder();
        
        typeBuilder.setName("Hachures");
        typeBuilder.setCRS(contours.getSchema().getCoordinateReferenceSystem());
        Class geomType = LineString.class;
        typeBuilder.add("Location", geomType);
        typeBuilder.add("Level", Double.class);
        
        final SimpleFeatureType TYPE = typeBuilder.buildFeatureType();
        GeometryFactory gFactory = JTSFactoryFinder.getGeometryFactory();
        SimpleFeatureBuilder fBuilder = new SimpleFeatureBuilder(TYPE);
        
        DefaultFeatureCollection collection = new DefaultFeatureCollection("internal",TYPE);
        
        int i = 0;
        for(CoordinateSequence hachure: hachuresDown){
            LineString line = gFactory.createLineString(hachure);
            fBuilder.add(line);
            fBuilder.add(hachureLevels.get(i));
            SimpleFeature feature = fBuilder.buildFeature(null);
            collection.add(feature);
            i++;
        }
        
        for(CoordinateSequence hachure: hachuresUp){
            LineString line = gFactory.createLineString(hachure);
            fBuilder.add(line);
            SimpleFeature feature = fBuilder.buildFeature(null);
            collection.add(feature);
        }
        
        return collection;
        
    }
    
    /**
     * Returns seed points that are used for hachure construction
     * @return 
     */
    public FeatureCollection getSeeds(){
        if(seeds == null)
            return null;
        
        SimpleFeatureTypeBuilder typeBuilder = new SimpleFeatureTypeBuilder();
        
        typeBuilder.setName("Seeds");
        typeBuilder.setCRS(contours.getSchema().getCoordinateReferenceSystem());
        Class geomType = Point.class;
        typeBuilder.add("Location", geomType);
        
        final SimpleFeatureType TYPE = typeBuilder.buildFeatureType();
        GeometryFactory gFactory = JTSFactoryFinder.getGeometryFactory();
        SimpleFeatureBuilder fBuilder = new SimpleFeatureBuilder(TYPE);
        
        DefaultFeatureCollection collection = new DefaultFeatureCollection();
        
        for(Coordinate c: seeds){
            Point point = gFactory.createPoint(c);
            fBuilder.add(point);
            SimpleFeature feature = fBuilder.buildFeature(null);
            collection.add(feature);
        }
        return collection;
    }
    
    /**
     * Traces hachures downslope and upslope, checks for intersections
     * @param zmin
     * @param zmax
     * @param nDown
     * @param nFirstPrev
     * @param nUp 
     */
    void getHachures(LineString cont, double curdist, double maxZ, 
                    double curZ, double minZ, int nFirstPrev){
        double L = 0; // line length
        double l ; // line segment length
        int nFirst = hachuresDown.size();
        int nFirstUp = hachuresUp.size();
        
        L = cont.getLength();
        
        // round number of steps
        
        double N = L/curdist;
        
        boolean toUp = (Math.ceil(N)-N) < (N-Math.floor(N));
        
        double newdist = (toUp) ? L/Math.ceil(N) : L/Math.floor(N);
        
        N = L/newdist;
        
        if(N<3) // if the contour is too short
            return;
        
        double s = newdist; // accumulated distance sum
        
        Coordinate prev = cont.getPointN(0).getCoordinate();
        Coordinate cur;
        
        for(int i = 1; i<cont.getNumPoints(); i++){
            cur = cont.getPointN(i).getCoordinate();
            l = AutolabMath.length(cur, prev);
            if(l==0)
                continue;
            s+=l;
            
            while(s>=newdist){
                double res = s - newdist; // residual
                
                l = AutolabMath.length(cur, prev);
                
                double t = 1 - res/l;
                
                if(t<0)
                    t = 0;
                
                Coordinate seedcoord = new Coordinate(
                        prev.x + (cur.x - prev.x)*t,
                        prev.y + (cur.y - prev.y)*t
                );
                
                if(!isVector){
                    // get downslope hachures
                    CoordinateSequence hachDown = getDownHachureFrom(seedcoord, minZ, nFirst,curdist);
                    if(hachDown != null)
                        if(hachDown.size()>1){
                           hachuresDown.add(hachDown);
                           hachureLevels.add(curZ);
                        }

                    // get upslope hachures
                    CoordinateSequence hachUp = getUpHachureFrom(seedcoord, maxZ, nFirst, nFirstPrev, nFirstUp, curdist);
                    if(hachUp != null)
                        if(hachUp.size()>1){
                           hachuresUp.add(hachUp);
                           hachureLevels.add(maxZ);
                        }
                } else {
                    // get stream hachures
                    CoordinateSequence hachDown = getStreamHachureFrom(seedcoord, minZ, curZ, maxZ, nFirst, curdist);
                    if(hachDown != null)
                        if(hachDown.size()>1){
                           hachuresDown.add(hachDown);
                           hachureLevels.add(curZ);
                        }
                }
                
                seeds.add(seedcoord);
                
                prev = seedcoord;
                s = res;
            }
            prev = cur;
        }
        return;   
    }
    
    /**
     * Iterates over each pair of adjacent hachures and inserts 
     * a new between them
     * @param zMin
     * @param n1
     * @param nFirst 
     */
    void insertHachures(double curdist, double zMin, double zCur, double zMax, int nFirstCont, int nFirst){
        int nLastCont = hachuresDown.size();
        if(!isVector)
            for(int i = nFirstCont+1; i < nLastCont-1; i++)
                insertDownHachures(curdist,i, i+1, 0, 0, zMin, 1, nFirst);
        else
            for(int i = nFirstCont+1; i < nLastCont-1; i++)
                insertStreamHachures(curdist, i, i+1, 0, 0, zMin, zCur, zMax, 1, nFirst);
    }
    
    /**
     * Inserts new hachures between n1 and n2 recursively 
     * until allowed depth is reached
     * @param n1 left hachure
     * @param n2 right hachure
     * @param k1 anchor vertex at left hachure
     * @param k2 anchor vertex at right hachure
     * @param zMin minimum z value at which the tracing is terminated
     * @param zMax maximum z value at which the tracing is terminated
     * @param depth current depth
     * @param nFirst first hachure number at current level
     */
    void insertStreamHachures(double curdist, int n1, int n2, 
                         int k1, int k2, 
                         double zMin, double zCur, double zMax, int depth, int nFirst){
        if(depth == maxDepth)
            return;
        
        Coordinate p1,p2,p3;
        
        // find how much points we can analyse
        int n = Math.min(hachuresDown.get(n1).size()-k1, 
                         hachuresDown.get(n2).size()-k2);
        
        for(int i = 0; i < n; i++){
            p1 = hachuresDown.get(n1).getCoordinate(i+k1);
            p2 = hachuresDown.get(n2).getCoordinate(i+k2);
            
            double dist = AutolabMath.length(p1, p2);
            
            if(dist > 2*curdist){ // if the space is too wide, then insert
                p3 = new Coordinate(0.5*(p1.x+p2.x),
                                    0.5*(p1.y+p2.y)); // in the middle
                CoordinateSequence hach = getStreamHachureFrom(p3,zMin,zCur,zMax,nFirst,curdist);
                if(hach != null){ // if the hachure has been created
                    hachuresDown.add(hach);
                    hachureLevels.add(zCur);
                    int nlast = hachuresDown.size()-1;
                    insertStreamHachures(curdist,n1,nlast,i+k1,0,zMin,zCur,zMax,depth+1,nFirst);
                    insertStreamHachures(curdist,n2,nlast,i+k2,0,zMin,zCur,zMax,depth+1,nFirst);
                    return;
                }
            }
        }
    }
    
    /**
     * Inserts new hachures between n1 and n2 recursively 
     * until allowed depth is reached
     * @param n1 left hachure
     * @param n2 right hachure
     * @param k1 anchor vertex at left hachure
     * @param k2 anchor vertex at right hachure
     * @param zMin minimum z value at which the tracing is terminated
     * @param depth current depth
     * @param nFirst first hachure number at current level
     */
    void insertDownHachures(double curdist, int n1, int n2, 
                            int k1, int k2, 
                            double zMin, int depth, int nFirst){
        if(depth == maxDepth)
            return;
        
        Coordinate p1,p2,p3;
        
        // find how much points we can analyse
        int n = Math.min(hachuresDown.get(n1).size()-k1, 
                         hachuresDown.get(n2).size()-k2);
        
        for(int i = 0; i < n; i++){
            p1 = hachuresDown.get(n1).getCoordinate(i+k1);
            p2 = hachuresDown.get(n2).getCoordinate(i+k2);
            
            double dist = AutolabMath.length(p1, p2);
            
            if(dist > 2*this.dist){ // if the space is too wide, then insert
                p3 = new Coordinate(0.5*(p1.x+p2.x),
                                    0.5*(p1.y+p2.y)); // in the middle
                CoordinateSequence hach = getDownHachureFrom(p3,zMin,nFirst,curdist);
                if(hach != null){ // if the hachure has been created
                    hachuresDown.add(hach);
//                    hachureLevels.add(zCur);
                    int nlast = hachuresDown.size()-1;
                    insertDownHachures(curdist,n1,nlast,i+k1,0,zMin,depth+1,nFirst);
                    insertDownHachures(curdist,n2,nlast,i+k2,0,zMin,depth+1,nFirst);
                    return;
                }
            }
        }
    }
    
    /**
     * Traces stream line from the point p forward until reaching 
     * lower or higher speed value
     * @param p start point
     * @param zMin minimum speed value
     * @param zMax maximum speed value
     * @param nDown
     * @return 
     */
    CoordinateSequence getStreamHachureFrom(Coordinate p, double zMin, double zlevel, double zMax, int nDown, double curdist){
        double[] G = new double[2]; // gradient vector
        Coordinate p1, p2, p3; // consecutive hachure vertices
        double zCur, zPre; // current and previous heights
        double t, l, A;
	int k;
        
        p1 = new Coordinate(p);
        
        int ij[]; // define the cell
        
        zPre = gp.getZxy(p.x, p.y);
        
        // outside model or lower that zMin
        if (zPre == Double.NaN || zPre < zMin) 
            return null;
        
        
        ArrayList<Coordinate> hach =  new ArrayList<>();
        
        hach.add(p1);
        
        // controlling circular streams!
        int startmonitoring = (int)Math.ceil(dist/step)*2;
        
        while(true){
            
            G[0] = gpx.getZxy(p1.x, p1.y);
            G[1] = gpy.getZxy(p1.x, p1.y);
            
            if (G[0] == Double.NaN || G[1] == Double.NaN)
                break;
            
            l = Math.sqrt(G[0]*G[0]+G[1]*G[1]);
            
            ij = gp.getIJfromXY(p1.x, p1.y);
            
            p2 = new Coordinate();
            
            p2.x = p1.x - step*G[0]/l;
            p2.y = p1.y - step*G[1]/l;
            
            if(hach.size()>startmonitoring && AutolabMath.length(p, p2) <= dist)
                break;
            
            zCur = gp.getZxy(p2.x, p2.y);
            
            // outside model or speed is too low
            if(zCur == Double.NaN || zCur < minSlope) 
                break;
            
            zPre = gp.getZxy(p1.x, p1.y);
            
            if(zCur < zMin) {
                t = (zCur - zMin) / (zCur - zPre);
                p3 = new Coordinate();
                p3.x = (1 - t)*p1.x + t*p2.x;
                p3.y = (1 - t)*p1.y + t*p2.y;
                hach.add(p2);
                break;
            } else if(zCur > zMax) {
                t = (zCur - zMax) / (zCur - zPre);
                p3 = new Coordinate();
                p3.x = (1 - t)*p1.x + t*p2.x;
                p3.y = (1 - t)*p1.y + t*p2.y;
                hach.add(p2);
                break;
            } else if((hach.size()>3) && ((zPre < zlevel)&&(zCur > zlevel))) {
                t = (zCur - zlevel) / (zCur - zPre);
                p3 = new Coordinate();
                p3.x = (1 - t)*p1.x + t*p2.x;
                p3.y = (1 - t)*p1.y + t*p2.y;
                hach.add(p2);
                break;
            } else if((hach.size()>3) && ((zPre > zlevel)&&(zCur < zlevel))) {
                t = (zCur - zlevel) / (zCur - zPre);
                p3 = new Coordinate();
                p3.x = (1 - t)*p1.x + t*p2.x;
                p3.y = (1 - t)*p1.y + t*p2.y;
                hach.add(p2);
                break;
            } else {
                p3 = new Coordinate(p2);
                p1 = new Coordinate(p2);
                hach.add(p2);
            }
            
            // check angle
            k = hach.size();
            if (k > 3)
                if(!checkAngle(hach.get(k-1),
                                hach.get(k-2),
                                hach.get(k-3))){
                    hach.remove(k-1);
                    break;
                }
        } // while
        
        // final check angle
        k = hach.size();
        if (k > 3)
            if(!checkAngle(hach.get(k-1),
                           hach.get(k-2),
                           hach.get(k-3))){
                hach.remove(k-1);
                k--;
            }
        
        if(k > 1){
            int n = checkDistances(hach,nDown,hachuresDown.size(),true, curdist);
            if(n > 0)
                for(; k > n; k--)
                    hach.remove(k-1);
        }
        
        
        Coordinate[] result = new Coordinate[hach.size()];
        
        if(k > 1)
            return new CoordinateArraySequence(hach.toArray(result));
        else 
            return null;
    }
    
    /**
     * Traces flow line from the point p downslope until reaching Zext value
     * @param p
     * @param zMin
     * @param nDown
     * @param nDown2
     * @param nUp
     * @return 
     */
    CoordinateSequence getDownHachureFrom(Coordinate p, double zMin, int nDown, double curdist){
        double[] G = new double[2]; // gradient vector
        Coordinate p1, p2, p3; // consecutive hachure vertices
        double zCur, zPre; // current and previous heights
        double t, l, A;
	int k;
        
        p1 = new Coordinate(p);
        
        int ij[]; // define the cell
        
        zPre = gp.getZxy(p.x, p.y);
        
        // outside model or lower that zMin
        if (zPre == Double.NaN || zPre < zMin) 
            return null;
        
        
        ArrayList<Coordinate> hach =  new ArrayList<>();
        
        hach.add(p1);
        
        while(true){
            
            G = gp.getGradXY(p1.x, p1.y);
            if (G == null)
                break;
            l = Math.sqrt(G[0]*G[0]+G[1]*G[1]);
            ij = gp.getIJfromXY(p1.x, p1.y);
            
            double slope = gp.getSlope(ij[0], ij[1]);
            
            if(slope <= minSlope)
                break;
            
            p2 = new Coordinate();
            
            p2.x = p1.x - step*G[0]/l;
            p2.y = p1.y - step*G[1]/l;
            
            zCur = gp.getZxy(p2.x, p2.y);
            
            if(zCur == Double.NaN) // outside model
                break;
            
            zPre = gp.getZxy(p1.x, p1.y);
            
            if(zCur < zMin){
                t = (zMin - zCur) / (zCur - zPre);
                p3 = new Coordinate();
                p3.x = (1 - t)*p1.x + t*p2.x;
                p3.y = (1 - t)*p1.y + t*p2.y;
                hach.add(p3);
                break;
            } else if(zCur < zPre){
                p3 = new Coordinate(p2);
                p1 = new Coordinate(p2);
                hach.add(p2);
            } else break;
            
            // check angle
            k = hach.size();
            if (k > 3)
                if(!checkAngle(hach.get(k-1),
                                hach.get(k-2),
                                hach.get(k-3))){
                    hach.remove(k-1);
                    break;
                }
        } // while
        
        // final check angle
        k = hach.size();
        if (k > 3)
            if(!checkAngle(hach.get(k-1),
                           hach.get(k-2),
                           hach.get(k-3))){
                hach.remove(k-1);
                k--;
            }
        
        if(k > 1){
            int n = checkDistances(hach,nDown,hachuresDown.size(),true,curdist);
            if(n > 0)
                for(; k > n; k--)
                    hach.remove(k-1);
        }
        
        
        Coordinate[] result = new Coordinate[hach.size()];
        
        if(k > 1)
            return new CoordinateArraySequence(hach.toArray(result));
        else 
            return null;
    }
    
    /**
     * Traces flow line from the point p upslope until reaching zMax value
     * @param p
     * @param Zext
     * @param nFirst
     * @param nFirstPrev
     * @param nFirstUp
     * @return 
     */
    CoordinateSequence getUpHachureFrom(Coordinate p, double zMax, int nFirst,
                               int nFirstPrev, int nFirstUp, double curdist){
        double[] G = new double[2]; // gradient vector
        Coordinate p1, p2, p3; // consecutive hachure vertices
        double zCur, zPre; // current and previous heights
        double t, l, A;
	int k;
        
        p1 = new Coordinate(p);
        
        int ij[]; // define the cell
        
        zPre = gp.getZxy(p.x, p.y);
        
        // outside model or higher that maximum
        if (zPre == Double.NaN || zPre > zMax) 
            return null;
        
        
        ArrayList<Coordinate> hach =  new ArrayList<>();
        
        hach.add(p1);
        
        while(true){
            
            G = gp.getGradXY(p1.x, p1.y);
            if (G == null)
                break;
            
            l = Math.sqrt(G[0]*G[0]+G[1]*G[1]);
            ij = gp.getIJfromXY(p1.x, p1.y);
            
            double slope = gp.getSlope(ij[0], ij[1]);
            
            if(slope <= minSlope)
                break;
            
            p2 = new Coordinate();
            
            p2.x = p1.x + step*G[0]/l;
            p2.y = p1.y + step*G[1]/l;
            
            zCur = gp.getZxy(p2.x, p2.y);
            
            if(zCur == Double.NaN) // outside model
                break;
            
            zPre = gp.getZxy(p1.x, p1.y);
            
            if(zCur >= zMax){
                hach.removeAll(hach);
                return null;
//                t = (zMax - zCur) / (zCur - zPre);
//                p3 = new Coordinate();
//                p3.x = (1 - t)*p1.x + t*p2.x;
//                p3.y = (1 - t)*p1.y + t*p2.y;
//                hach.add(p3);
//                break;
            } else if(zCur > zPre){
                p3 = new Coordinate(p2);
                p1 = new Coordinate(p2);
                hach.add(p2);
            } else break;
            
            // check angle
            k = hach.size();
            if (k > 3)
                if(!checkAngle(hach.get(k-1),
                                hach.get(k-2),
                                hach.get(k-3))){
                    hach.remove(k-1);
                    break;
                }
        } // while
        
        // final check angle
        k = hach.size();
        if (k > 3)
            if(!checkAngle(hach.get(k-1),
                           hach.get(k-2),
                           hach.get(k-3))){
                hach.remove(k-1);
                k--;
            }
        
        if(k > 1){
            int n = checkDistances(hach,nFirstUp,hachuresUp.size(),false, curdist);
            if(n > 0)
                for(; k > n; k--)
                    hach.remove(k-1);
            
            n = checkDistances(hach,nFirstPrev,nFirst,true, curdist);
            if(n > 0)
                for(; k > n; k--)
                    hach.remove(k-1);
        }
        
        Coordinate[] result = new Coordinate[hach.size()];
        
        if(k > 1)
            return new CoordinateArraySequence(hach.toArray(result));
        else 
            return null;
    }
    
    
    /**
     * Checks an angle between three consecutive line points
     * @param p1
     * @param p2
     * @param p3
     * @return 
     */
    boolean checkAngle(Coordinate p1,
                       Coordinate p2,
                       Coordinate p3){
        double num = (p2.x - p1.x)*(p3.x - p2.x) + (p2.y - p1.y)*(p3.y - p2.y);
	double den = Math.pow((p2.x - p1.x)*(p2.x - p1.x) + 
                              (p2.y - p1.y)*(p2.y - p1.y),0.5)*
                     Math.pow((p3.x - p2.x)*(p3.x - p2.x) + 
                              (p3.y - p2.y)*(p3.y - p2.y),0.5);
	if(den>0.0){
		double cosA = num/den;
		return (cosA >= Math.cos(maxTurn));
	}
        
        return false;
    }
    
    /**
     * Checks distances and intersections between current (hach) hachure and 
     * those previously created;
     * @param hach
     * @param n
     * @param down
     * @return 
     */
    int checkDistances(ArrayList<Coordinate> hach, int n1, int n2, boolean down, double curdist){
        double length;
        boolean intersect;
        
        ArrayList<CoordinateSequence> hachures = (down) ? hachuresDown : hachuresUp;
        
        for(int i = 1; i < hach.size(); i++){
            for(int j = n1; j < n2; j++){
                CoordinateSequence hachref = hachures.get(j);
                int npts = hachref.size();
                for(int k = 1; k < npts; k++){
                    length = AutolabMath.length(hach.get(i), hachref.getCoordinate(k));
                    if(length < 0.5*curdist)
                        return i;
                    intersect = AutolabMath.testIntersection(
                                            hach.get(i-1),
                                            hach.get(i), 
                                            hachref.getCoordinate(k-1), 
                                            hachref.getCoordinate(k));
                    if(intersect)
                        return i;
                }
            }
        }
        return 0;
    }
}
