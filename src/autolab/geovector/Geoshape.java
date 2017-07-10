/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.geovector;
import autolab.symbols.Symbol;
import java.awt.*;
import java.util.ArrayList;
import javafx.util.Pair;

/**
 *
 * @author timofey
 */
public abstract class Geoshape {
    String label;
    double minX, maxX;
    double minY, maxY;
    double centerX, centerY;
    ArrayList<Pair> attributes;
    
    public abstract void paint(Graphics2D g2d, Symbol symbol);
    protected abstract void updateBounds();
    public Geoshape(){
        label = new String("Yaya");
        minX = minY = Double.MAX_VALUE;
        maxX = maxY = Double.MIN_VALUE;
    }
    public void setLabel(String s){
        label = s;
    }
    public String getLabel(){
        return label;
    }
    public double getCenterX(){
        return centerX;
    }
    public double getCenterY(){
        return centerX;
    }
    public double getMinX(){
        return minX;
    }
    public double getMaxX(){
        return maxX;
    }
    public double getMinY(){
        return minY;
    }
    public double getMaxY(){
        return maxY;
    }
    public void updateCenter(){
        centerX = 0.5*(minX + maxX);
        centerY = 0.5*(minY + maxY);
    }
}
