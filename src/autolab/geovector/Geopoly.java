/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.geovector;
import autolab.symbols.PolySimple;
import autolab.symbols.Symbol;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.awt.*;

/**
 *
 * @author timofey
 */
public class Geopoly extends Geoshape {
    private ArrayList<Double> X;
    private ArrayList<Double> Y;
    
    public Geopoly(){
        super();
        X = new ArrayList<>();
        Y = new ArrayList<>();
    }
    
    public Geopoly(double[] Xcoords, double[] Ycoords){
        super();
        X = new ArrayList(Xcoords.length);
        Y = new ArrayList(Ycoords.length);
        for (int i = 0; i<X.size(); i++){
            X.set(i,Xcoords[i]);
            Y.set(i,Ycoords[i]);   
        }
        updateBounds();
    }
    
    @Override
    protected void updateBounds(){
        Double x,y;
        for (int i = 0; i < X.size(); i++){
            x = X.get(i);
            y = Y.get(i);
            if (x < minX){
                minX = x;
            }
            else if (x > maxX){
                maxX = x;
            }
            if (y < minY){
                minY = y;
            }
            else if (y > maxY){
                maxY = y;
            }
        }
        updateCenter();
    }
    
    public void addPoint(Double x, Double y){
        X.add(x);
        Y.add(y);
        updateBounds();
    }
    
    public void delPoint(int i){
        X.remove(i);
        Y.remove(i);
        updateBounds();
    }
    
    public void insPoint(int i, Double x, Double y){
        X.add(i, x);
        Y.add(i, y);
        updateBounds();
    }
    
    public void setPoint(int i, Double x, Double y){
        X.set(i, x);
        Y.set(i, y);
        updateBounds();
    }
    
    public void setSize(int n){
        X.ensureCapacity(n);
        Y.ensureCapacity(n);
    }
    public int[] getXarray(){
        int xarr[] = new int[X.size()];
        for (int i = 0; i < X.size(); i++){
            xarr[i] = (int)(double)X.get(i);
        }
        return xarr;
    }
    
    public int[] getYarray(){
        int yarr[] = new int[Y.size()];
        for (int i = 0; i < Y.size(); i++){
            yarr[i] = (int)(double)Y.get(i);
        }
        return yarr;
    }
    
    @Override
    public void paint(Graphics2D g2d, Symbol symbol){
        PolySimple s = (PolySimple)symbol;
        int xarr[] = getXarray();
        int yarr[] = getYarray();
        g2d.setColor(s.getFillColor());
        g2d.fillPolygon(xarr, yarr, xarr.length);
        for (int i = 0; i < s.getStrNumber(); i++){
             g2d.setColor(s.getStrokeColor(i));
             g2d.setStroke(s.getStroke(i));
             g2d.drawPolygon(xarr, yarr, xarr.length);
        }
    }
    public void paintFill(Graphics2D g2d, Color fill){
        int xarr[] = getXarray();
        int yarr[] = getYarray();
        g2d.setColor(fill);
        g2d.fillPolygon(xarr, yarr, xarr.length);
    };
    
    public void paintStroke(Graphics2D g2d, BasicStroke stroke, Color color){
        int xarr[] = getXarray();
        int yarr[] = getYarray();
        g2d.setColor(color);
        g2d.setStroke(stroke);
        g2d.drawPolygon(xarr, yarr, xarr.length);
    };
    
    public void paintLabel(Graphics2D g2d, Color color){
        g2d.setColor(color);
        g2d.drawString(label, (float)centerX, (float)centerY);
    };
}
