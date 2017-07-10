/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.geovector;
import java.awt.Graphics2D;
import java.util.ArrayList;
import autolab.symbols.Symbol;
/**
 *
 * @author Timofey Samsonov, Department of Cartography and Geoinformatics
 * Lomonosov Moscow State University, Russia.
 */
public abstract class Geoset {
    /**
     * Should be implemented by inherited classes
     * @param g2d Graphics context
     */
    Symbol symbol;
    ArrayList<Geoshape> set;
    boolean drawLabels;
    
    public abstract void paint(Graphics2D g2d);
    
    public Geoset(){
        drawLabels = false;
    }
    
    public void addShape(Geoshape shape){
        set.add(shape);
    }
    public void setJoinStrokes(boolean isJoined){
        symbol.setJoinStrokes(isJoined);
    }
    public void setSymbol(Symbol s){
        symbol = s; 
    };
    public Symbol getSymbol(){
        return symbol;
    }
    public void setLabel(int index, String s){
        set.get(index).setLabel(s);
    }
    public String getLabel(int index){
        return set.get(index).getLabel();
    }
    public boolean getDrawLabels(){
        return drawLabels;
    }
    public void setDrawLabels(boolean draw){
        drawLabels = draw;
    }
}
