/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.geovector;

import java.awt.Graphics2D;
import java.util.*;
import autolab.symbols.Symbol;

/**
 *
 * @author timofey
 */
public class Geoline extends Geoshape {
    ArrayList<Double> X;
    ArrayList<Double> Y;
    
    public Geoline(){
        X = new ArrayList<Double>();
        Y = new ArrayList<Double>();
    }
    
    public Geoline(double[] Xcoords, double[] Ycoords){
        X = new ArrayList(Xcoords.length);
        Y = new ArrayList(Ycoords.length);
        for (int i = 0; i<X.size(); i++){
            X.set(i,Xcoords[i]);
            Y.set(i,Ycoords[i]);   
        }
    }
    
    public Geopoint getPoint(int i){
        return new Geopoint(X.get(i),Y.get(i));
    }
    
    public void addPoint(Double x, Double y){
        X.add(x);
        Y.add(y);
    }
    
    public void addPoint(Geopoint p){
        X.add(p.x);
        Y.add(p.y);
    }
    
    public void delPoint(int i){
        X.remove(i);
        Y.remove(i);
    }
    
    public void insPoint(int i, Double x, Double y){
        X.add(i, x);
        Y.add(i, y);
    }
    
    public void setPoint(int i, Double x, Double y){
        X.set(i, x);
        Y.set(i, y);
    }
    
    public void setSize(int n){
        X.ensureCapacity(n);
        Y.ensureCapacity(n);
    }
    
    public int getSize(){
        return X.size();
    }
    
    public double getLength(){
        return 0;
    }

    @Override
    public void paint(Graphics2D g2d, Symbol symbol) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected void updateBounds() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    
}
