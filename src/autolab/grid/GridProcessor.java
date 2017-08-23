/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.grid;

import java.awt.image.renderable.RenderableImage;
import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;

/**
 *
 * @author tsamsonov
 */
public class GridProcessor {
    
    private double Z[][];
    private GridHeader h;
    
    
    public GridProcessor(GridCoverage2D gc){
        h = new GridHeader(gc);
        
        
        System.out.println("Converting raster");
        System.out.println(h.nrow);
        System.out.println(h.ncol);
        Z = new double[h.ncol][h.nrow];
        for(int j = 0; j < h.nrow; j++){
            for(int i = 0; i < h.ncol; i++){
                GridCoordinates2D coords = new GridCoordinates2D(i,j);
                double values[] = new double[1];
                Z[i][h.nrow-j-1] = gc.evaluate(coords, values)[0];
            }
        }
    }
    
    public GridProcessor(){
        
    }
    
    public GridHeader getHeader(){
        return h;
    }
    
    /**
     * Function returns XY coordinates of the (i,j) cell
     * @param header DEM header
     * @param i the row
     * @param j the column
     * @return 
     */
    public double[] getXYfromIJ(int i, int j){
        double coords[] = new double[2];
        
        if(i<0 || i > h.ncol || j<0 || j > h.nrow){
            return null;
        }
        
        coords[0] = i*h.res + h.xmin;
        coords[1] = j*h.res + h.ymin;
        
        return coords;
    }
    
    /**
     * Returns (i,j) of the cell in which (x,y) coordinates fall to
     * @param h
     * @param x
     * @param y
     * @return 
     */
    public int[] getIJfromXY(double x, double y){
        int coords[] = new int[2];
        
        if(x<=h.xmax && x >= h.xmin && y <= h.ymax && y >= h.ymin){
            double dx = x-h.xmin;
            double dy = y-h.ymin;
            
            if (x==h.xmax){
                coords[0] = h.ncol-2;
            } else {
                coords[0] = (int) Math.floor(dx/h.res);
            }
            
            if (y==h.ymax){
                coords[1] = h.nrow-2;
            } else {
                coords[1] = (int) Math.floor(dy/h.res);
            }
            return coords;
        }
        else return null;
    }
    
    /**
     * Function returns Z value interpolated at (x,y) using bilinear function
     * @param Z matrix
     * @param x coordinate
     * @param y coordinate
     * @return 
     */
    public double getZxy(double x, double y){
        int i,j;
        double x1,y1,x2,y2;
        int[] ij = getIJfromXY(x,y);
        if(ij == null)
            return Double.NaN;
        i=ij[0];
        j=ij[1];
        if(i == h.ncol-1 || j == h.nrow-1)
            return Double.NaN;
        if (ij != null){
            double[] xy1 = getXYfromIJ(i,j);
            if(xy1 == null){
                return Double.NaN;
            }
            double dx = (x-xy1[0])/h.res;
            double dy = (y-xy1[1])/h.res;
            double Ax = Z[i+1][j] - Z[i][j];
            double Ay = Z[i][j+1] - Z[i][j];
            double Axy = Z[i][j] + Z[i+1][j+1] - Z[i+1][j] - Z[i][j+1];
            double Z0 = Z[i][j] + Ax*dx + Ay*dy + Axy*dx*dy;
            return Z0;
        } else {
            return Double.NaN;
        }
    }
    
    /**
     * Returns 3D normal to the surface at (i,j) cell
     * @param i coordinate
     * @param j coordinate
     * @return 
     */
    public double[] getNorm(int i, int j){
        double k = 1.0;
        double l = 1.0;
        double[] N = new double[3];
        
        if(i>0 && i<h.ncol-1){
            k = 2.0;
            N[0] = -2*h.zk*(Z[i+1][j] - Z[i-1][j])*h.res;
        } else if(i==0){
            N[0] = -h.zk*(Z[i+1][j] - Z[i][j])*h.res;
        } else {
            N[0] = -h.zk*(Z[i][j] - Z[i-1][j])*h.res;
        }
        
        if(j>0 && j < h.nrow-1){
            l = 2.0;
            N[1] = -2*h.zk*(Z[i][j+1] - Z[i][j-1])*h.res;
        } else if(j==0){
            N[1] = -h.zk*(Z[i][j+1] - Z[i][j])*h.res;
        } else {
            N[1] = -h.zk*(Z[i][j] - Z[i][j-1])*h.res;
        }
        
        N[2] = k*l*h.res*h.res;
        return N;
    }
    
    /**
     * Returns 2D gradient of the surface at (i,j) cell
     * @param i
     * @param j
     * @param Z
     * @param h
     * @return 
     */
    public double[] getGrad(int i, int j){
        double[] G = new double[2];
        
        if(i>0 && i<h.ncol-1){
            G[0] = -0.5*h.zk*(Z[i+1][j] - Z[i-1][j])/h.res;
        } else if(i==0){
            G[0] = -h.zk*(Z[i+1][j] - Z[i][j])/h.res;
        } else {
            G[0] = -h.zk*(Z[i][j] - Z[i-1][j])/h.res;
        }
        
        if(j>0 && j < h.nrow-1){
            G[1] = -0.5*h.zk*(Z[i][j+1] - Z[i][j-1])/h.res;
        } else if(j==0){
            G[1] = -h.zk*(Z[i][j+1] - Z[i][j])/h.res;
        } else {
            G[1] = -h.zk*(Z[i][j] - Z[i][j-1])/h.res;
        }
        
        return G;
    }
    
    /**
     * Calculates gradient from the point (x,y)
     * @param x
     * @param y
     * @return 
     */
    public double[] getGradXY(double x, double y){
        int i,j;
        double x1,y1,x2,y2;
        int[] ij = getIJfromXY(x,y);
        if(ij == null)
            return null;
        i=ij[0];
        j=ij[1];
        
        if(i == h.ncol-1 || j == h.nrow-1)
            return null;
        
        if (ij != null){
            double[] xy1 = getXYfromIJ(i,j);
            if(xy1 == null){
                return null;
            }
            double dx = (x-xy1[0])/h.res;
            double dy = (y-xy1[1])/h.res;
            double Ax = Z[i+1][j] - Z[i][j];
            double Ay = Z[i][j+1] - Z[i][j];
            double Axy = Z[i][j] + Z[i+1][j+1] - Z[i+1][j] - Z[i][j+1];
            double[] G = new double[2];
            G[0] = h.zk*(Ax + dy*Axy);
            G[1] = h.zk*(Ay + dx*Axy);
            return G;
        } else {
            return null;
        }
    }
    
    /**
     * Returns surface slope in (i,j) cell
     * @param i
     * @param j
     * @param Z
     * @param h
     * @return 
     */
    public double getSlope(int i, int j){
        double[] N = new double[3];
        N = getNorm(i,j);
        if (N != null){
            double len = Math.pow(N[0]*N[0] + N[1]*N[1] + N[2]*N[2], 0.5);
            double slope = Math.acos(N[2]/len);
            return slope;
        } else {
            return Double.NaN;
        }
    } 
}
