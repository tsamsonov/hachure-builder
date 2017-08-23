/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.math;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Point;

/**
 *
 * @author tsamsonov
 */
public class AutolabMath {
    /**
     * Calculates logarithm by base n.
     * @param x — value for which logarithm should be calculated
     * @param n — base of the logarithm
     * @return 
     */
    public static double logn(double x, double n){
        return java.lang.Math.log(x)/java.lang.Math.log(n);
    }
    
    public static double azimuth(Coordinate p1, Coordinate p2){
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        
        double a;
        
        if (dx != 0 && dy != 0){
            a = java.lang.Math.atan(dy/dx);
            if (dx < 0 && dy > 0){ // second quarter
                a = -a + java.lang.Math.PI*0.5;
            } else if (dx > 0){    // third and fourth quarter
                a = java.lang.Math.PI*1.5 + a;
            }
        } else if (dx == 0 && dy > 0){
            a = java.lang.Math.PI;
        } else if (dx == 0 && dy < 0){
            a = 0;
        } else {
            return -1;
        }
        
        return a;
    }
    
    // mathematic azimuth
    public static double azimuth(double dx, double dy){
        double a;
        
        if (dx != 0 && dy != 0){
            a = Math.atan(dy/dx);
            if (dx < 0){ // second & third quarter
                a = Math.PI + a;
            } else if (dy < 0 && dx > 0){    // third and fourth quarter
                a = Math.PI*2 + a;
            }
        } else if (dx == 0 && dy > 0){
            a = Math.PI;
        } else if (dx == 0 && dy < 0){
            a = 0;
        } else if (dx > 0 && dy == 0){
            a = Math.PI * 0.5;
        } else if (dx < 0 && dy == 0){
            a = Math.PI * 1.5;
        } else return -1;
        
        return a;
    }
    
    public static double toGeographicAzimuth(double angle){
        double geographic = 0.0;
        
        if (angle >=0 && angle <= Math.PI*0.5){
            geographic = 0.5*Math.PI - angle;
        } else {
            geographic = 2.5*Math.PI - angle;
        }
        
        return geographic;
    }
    
    public static double toMathematicAzimuth(double geographic){
        double angle = 0.0;
        
        if (geographic >=0 && geographic <= Math.PI*0.5){
            angle = 0.5*Math.PI - geographic;
        } else {
            angle = 2.5*Math.PI - geographic;
        }
        
        return angle;
    }
    
    public static double averageAzimuth(double a1, double a2){
        double avg = (a1+a2)/2;
        double delta = Math.abs(a1 - a2);
        if (delta > 180){
            avg = avg + 180;
            if (avg > 360)
                avg -= 360;
        }
        return avg;
    }
    
    /**
     * Calculates the length of the line segment between points p1 and p2
     * @param p1
     * @param p2
     * @return 
     */
    public static double length(Point p1, Point p2){
        return Math.sqrt((p1.getX()-p2.getX())*(p1.getX()-p2.getX()) + 
                         (p1.getY()-p2.getY())*(p1.getY()-p2.getY()));
    }
    /**
     * Calculates the length of the line segment between points p1 and p2
     * @param c1
     * @param c2
     * @return 
     */
    public static double length(Coordinate c1, Coordinate c2){
        return Math.sqrt((c1.x-c2.x)*(c1.x-c2.x) + 
                         (c1.y-c2.y)*(c1.y-c2.y));
    }
    
    /**
     * Defines from which side of p1-p2 line segment p3 point is located
     * @param p1
     * @param p2
     * @param p3
     * @return 
     */
    public static boolean pointSide(Coordinate p1, Coordinate p2, Coordinate p3){
        double L = (p3.x - p1.x)*(p2.y - p1.y) - (p3.y - p1.y)*(p2.x - p1.x);
        return L<0;
    }
    
    
    /**
     * Tests intersection between p1-p2 and p3-p4 line segments
     * @param p1
     * @param p2
     * @param p3
     * @param p4
     * @return 
     */
    public static boolean testIntersection(Coordinate p1, Coordinate p2, 
                                           Coordinate p3, Coordinate p4){
        if((pointSide(p1,p2,p3)!=pointSide(p1,p2,p4)) &&
           (pointSide(p3,p4,p1)!=pointSide(p3,p4,p2)))
            return true;
        else return false;
    }
}
