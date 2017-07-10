/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.geovector;

import java.awt.Graphics2D;
import autolab.symbols.Symbol;

/**
 *
 * @author timofey
 */
public class Geopoint extends Geoshape {
    public double x;
    public double y;
    
    public Geopoint(){
        x = 0;
        y = 0;
    }
    
    public Geopoint(double x, double y){
        this.x = x;
        this.y = y;
    }
    
    // initializes from point p
    public Geopoint(Geopoint p){
        x = p.x;
        y = p.y;
    }
    
    // creates from two points p1 and p2
    public Geopoint(Geopoint p1, Geopoint p2){
        x = 0.5*(p1.x+p2.x);
        y = 0.5*(p1.y+p2.y);
    }
    
    void setX(double x, double y){
        this.x = x;
        this.y = y;
    }

    @Override
    public void paint(Graphics2D g2d, Symbol symbol) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected void updateBounds() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    // initializes from point p
    public void initFrom(Geopoint p){
        x = p.x;
        y = p.y;
    }
    
    // initializes from two points p1 and p2
    public void initFrom(Geopoint p1, Geopoint p2){
        x = 0.5*(p1.x+p2.x);
        y = 0.5*(p1.y+p2.y);
    }
}
