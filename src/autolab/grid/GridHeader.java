/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.grid;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;

/**
 *
 * @author tsamsonov
 */
public class GridHeader {
     
    public double xmin, xmax;
    public double ymin, ymax;
    public double zmin, zmax;
    public double res;
    public int nrow, ncol;
    public double noData = Float.NaN;
    public double zk;
    
    public GridHeader(){
        
    }
    /**
     * Grid header short constructor. Other values are calculated from parameters
     * @param x
     * @param y
     * @param rows
     * @param cols
     * @param cellsize 
     */
    public GridHeader(double x, double y, int rows, int cols, 
                        float z1, float z2, double cellsize){
        // Source parameters
        xmin = x;
        ymin = y;
        nrow = rows;
        ncol = cols;
        res = cellsize;
        zmin = z1;
        zmax = z2;
        
        // Calculated parameters
        xmax = xmin + cols*cellsize;
        ymax = ymin + rows*cellsize;
        
        // Default parameters
        noData = Float.NaN;
        zk = 1.0f;
    }
    
    public GridHeader(GridHeader h){
        // Source parameters
        xmin = h.xmin;
        ymin = h.ymin;
        nrow = h.nrow;
        ncol = h.ncol;
        res = h.res;
        zmin = h.zmin;
        zmax = h.zmax;
        
        // Calculated parameters
        xmax = xmin + h.ncol * res;
        ymax = ymin + h.nrow * res;
        
        // Default parameters
        noData = Float.NaN;
        zk = 1.0f;
    }
    
    public GridHeader(GridCoverage2D cov){
            zmin = (float) cov.getSampleDimension(0).getMinimumValue();
            zmax = (float) cov.getSampleDimension(0).getMaximumValue();
            zk = (float) cov.getSampleDimension(0).getScale();

            GridGeometry2D g2d = cov.getGridGeometry();
            nrow = g2d.getGridRange2D().height;
            ncol = g2d.getGridRange2D().width;

            res = g2d.getEnvelope2D().width/ncol;

            xmin = g2d.getEnvelope2D().getMinX()+res/2;
            xmax = g2d.getEnvelope2D().getMaxX()+res/2;

            ymin = g2d.getEnvelope2D().getMinY()+res/2;
            ymax = g2d.getEnvelope2D().getMaxY()+res/2;
    }
}
