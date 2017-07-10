/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.geovector;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.ArrayList;
import autolab.symbols.PolySimple;
import autolab.symbols.Symbol;

/**
 *
 * @author timofey
 */
public class Geopolyset extends Geoset {
    
     
    public Geopolyset(){
        set = new ArrayList<>();
        symbol = new PolySimple();
    }
    
    @Override
    public void setSymbol(Symbol s){
        this.symbol = (PolySimple)s;
    }
    
    @Override
    public void paint(Graphics2D g2d){
        if(!symbol.isJoinStrokes()) {
            for (Geoshape shp:set) {
                shp.paint(g2d, symbol);
            } 
        } 
        else{
            PolySimple ps = (PolySimple)symbol;
            for (int i = 0; i < ps.getStrNumber(); i++){
                for (Geoshape shp:set){
                    Geopoly poly  = (Geopoly)shp;
                    poly.paintStroke(g2d, ps.getStroke(i), ps.getStrokeColor(i));
                } 
            }
            for (Geoshape shp:set){
                Geopoly poly  = (Geopoly)shp;
                poly.paintFill(g2d, ps.getFillColor());
            }
        }
        if (drawLabels){
            for (Geoshape shp:set){
                Geopoly poly  = (Geopoly)shp;
                poly.paintLabel(g2d, Color.BLACK);
            }
        }       
    }
}
