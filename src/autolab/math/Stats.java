/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.math;

import java.util.ArrayList;
import java.lang.Math;
import java.util.Arrays;
import java.util.Collections;

/**
 *
 * @author tsamsonov
 */

class Pair<X, Y> implements Comparable{
    public final X x;
    public final Y y;
    
    public Pair(X x, Y y){
        this.x = x;
        this.y = y;
    }

    @Override
    public int compareTo(Object o) {
        if (x instanceof Double && o instanceof Double){
            Double d = (Double) o;
            Double i = (Double) x;
            if (i > d) 
                return 1;
            else if (i < d) 
                return -1;
            else 
                return 0;
        }
        return -2;
    }
    
}
public class Stats {
    
    // Methods used for reclassification of numerical arrays
    public enum ReclassMethod {
        DEFINED, 
        EQUAL, 
        QUANTILE, 
        STANDARD, 
        NATURAL, 
        GEOMETRIC, 
        NESTED, 
        HEADTAIL,
        EXPANDING;
        
        @Override
        public String toString(){
            switch(this){
                case DEFINED: return "Defined Intervals";
                case EQUAL: return "Equal Intervals";
                case QUANTILE: return "Quantile";
                case STANDARD: return "Standard Deviations";
                case NATURAL: return "Natural Breaks (Jenks)";
                case GEOMETRIC: return "Geometric (Smart Quantiles)";
                case NESTED: return "Nested Means";
                case HEADTAIL: return "Head/Tail Breaks";
                case EXPANDING: return "Expanding Quantiles";
                default: throw new IllegalArgumentException();
            }
        }
    }
    
    // Resampling methods
    
    public enum ResampleMethod{
        NEAREST,
        NATURAL,
        BILINEAR,
        BICUBUC,
        MAJORITY;
        
        @Override
        public String toString(){
            switch(this){
                case NEAREST: return "Nearest Neighbor";
                case NATURAL: return "Natural Neighbor";
                case BILINEAR: return "Bilinear";
                case BICUBUC: return "Bicubic";
                case MAJORITY: return "Majority";
                default: throw new IllegalArgumentException();
            }
        }
    }
    
    // Methods used for filtering
    public enum FilterMethod {
        VECTOR,
        SCALAR,
        MIN,
        MAX,
        MEAN,
        MEDIAN,
        MODE,
        STDEV,
        RANGE,
        GAUSSIAN,
        ROBERTS,
        SOBEL,
        PREWITT,
        CANNY,
        MANUAL;
        
        @Override
        public String toString(){
            switch(this){
                case VECTOR: return "Vector Average";
                case SCALAR: return "Scalar average";
                case MIN: return "Minimum";
                case MAX: return "Maximum";
                case MEAN: return "Mean";
                case MEDIAN: return "Median";
                case MODE: return "Mode";
                case STDEV: return "Standard deviations";
                case RANGE: return "Range";
                case GAUSSIAN: return "Gaussian";
                case ROBERTS: return "Roberts Edge Detection";
                case SOBEL: return "Sobel Edge Detection";
                case PREWITT: return "Prewitt Edge Detectin";
                case CANNY: return "Canny Edge Detection";
                case MANUAL: return "Manual";
                default: throw new IllegalArgumentException();
            }
        }
    }
    
     // Methods used for filtering
    public enum AggregateMethod {
        MIN,
        MAX,
        MEAN,
        MEDIAN,
        MODE,
        STDEV,
        RANGE;
        
        @Override
        public String toString(){
            switch(this){
                case MIN: return "Minimum";
                case MAX: return "Maximum";
                case MEAN: return "Mean";
                case MEDIAN: return "Median";
                case MODE: return "Mode";
                case STDEV: return "Standard deviations";
                case RANGE: return "Range";
                default: throw new IllegalArgumentException();
            }
        }
    }
    
    /**
     * Calculates statistical modes in values[] one-dimensional array. 
     * @param values array to be analysed
     * @param step step using for histogram calculation
     * @param quantity number of modes to be calculated. modes will be sorted by their frequencies 
     * @param ring true if codomain is closed. 
     * @return 
     */
    
    public static double[] Modes(double[] values, double step, int quantity, boolean ring){
        
        // copy values to new sorted array — needed for easy histogram calculation
        double valuesSorted[] = Arrays.copyOf(values, values.length);
        Arrays.sort(valuesSorted);
        
        double min = valuesSorted[0];
        double max = valuesSorted[values.length-1];
        
        // calculate the number of classes
        int nClasses = (int) ((max-min) / step) + (int) Math.ceil((max-min) % step);
        
        double freq[] = new double[nClasses]; // frequencies array
        
        
        
        // reset all frequencies to zero
        for (int i = 0; i < freq.length; i++){
            freq[i]=0.0;
        }
        
        // ceil is upper boundary for current class 
        double ceil = min + step;
                    
        int j = 0;

        // calculate frequencies
        for (int i = 0; i < freq.length; i++){
            if (j<values.length){
                while (values[j] < ceil){
                    freq[i]++;
                    j++;
                    if (j == values.length) break;
                }
            }
            ceil += step;
        }
        
        ArrayList<Pair> peaks = new ArrayList<>();
        
        // find all internal peaks
        for (int i = 1; i < freq.length-1; i++) {
            if ( (freq[i] > freq[i-1]) && (freq[i] > freq[i+1])){
                peaks.add(new Pair(freq[i], i));
            }
        }
        
        
        // find boundary peaks
        j = freq.length-1;
        if (ring){ // ring - for example, azimuth
            if ( (freq[j] > freq[j-1]) && (freq[j] > freq[0])){
                peaks.add(new Pair(freq[j], j));
            } else if ( (freq[0] > freq[j]) && (freq[0] > freq[1])){
                peaks.add(new Pair(freq[0], 0));
            }
        } else {
            if (freq[j] > freq[j-1]){
                peaks.add(new Pair(freq[j], j));
            }
            if (freq[0] > freq[1]){
                peaks.add(new Pair(freq[0], 0));
            }
        }
        
        // sort peaks by frequencies
        Collections.sort(peaks);
        
        /* if number of peaks is less than number of modes to be found, then 
        * function will return all the peaks, otherwise it will return
        * firs <quantity> peaks
        */
        int nModes = (peaks.size() < quantity) ? peaks.size() : quantity;
        
        double modes[] = new double[nModes]; // modes array
        
        // calculate mode values
        
        double m0, a;
        double b = 0, c = 0;
        
        for (int i = 0; i < nModes; i++){
            int modeIndex = (int)peaks.get(i).y;
            m0 = min + step*modeIndex;
            a = (double)peaks.get(i).x;
            
            // select b and c coefficients
            if ((modeIndex > 0) && (modeIndex < nClasses-1)){ // entry in between
                b = freq[modeIndex-1];
                c = freq[modeIndex+1];
            } else if (modeIndex == 0){ // first entry
                c = freq[modeIndex+1];
                b = (ring == true) ? freq[freq.length-1] : 0;
            } else if (modeIndex == freq.length-1) { // last entry
                b = freq[modeIndex-1];
                c = (ring == true) ? freq[0] : 0;
            }
            // calculate modes
            modes[i] = m0 + step*(a-b)/(2*a-b-c);
        }
        
        
        return modes;
    }
    
    public static double[] Modes(double[] values, int quantity, boolean ring){
        
        // copy values to new sorted array — needed for easy histogram calculation
        double valuesSorted[] = Arrays.copyOf(values, values.length);
        Arrays.sort(valuesSorted);
        
        double min = valuesSorted[0];
        double max = valuesSorted[values.length-1];
        
        double step = meanDelta(valuesSorted)*5;
        
        int freq[] = Histogram(values, step, false); // frequencies array
        
        int nClasses = freq.length;
        
        
        ArrayList<Pair> peaks = new ArrayList<>();
        
        // find all internal peaks
        for (int i = 1; i < freq.length-1; i++) {
            if ( (freq[i] > freq[i-1]) && (freq[i] > freq[i+1])){
                peaks.add(new Pair(freq[i], i));
            }
        }
        
        
        // find boundary peaks
        int j = freq.length-1;
        if (ring){ // ring - for example, azimuth
            if ( (freq[j] > freq[j-1]) && (freq[j] > freq[0])){
                peaks.add(new Pair(freq[j], j));
            } else if ( (freq[0] > freq[j]) && (freq[0] > freq[1])){
                peaks.add(new Pair(freq[0], 0));
            }
        } else {
            if (freq[j] > freq[j-1]){
                peaks.add(new Pair(freq[j], j));
            }
            if (freq[0] > freq[1]){
                peaks.add(new Pair(freq[0], 0));
            }
        }
        
        // sort peaks by frequencies
        Collections.sort(peaks);
        
        /* if number of peaks is less than number of modes to be found, then 
        * function will return all the peaks, otherwise it will return
        * firs <quantity> peaks
        */
        int nModes = (peaks.size() < quantity) ? peaks.size() : quantity;
        
        double modes[] = new double[nModes]; // modes array
        
        // calculate mode values
        
        double m0, a;
        double b = 0, c = 0;
        
        for (int i = 0; i < nModes; i++){
            int modeIndex = (int)peaks.get(i).y;
            m0 = min + step*modeIndex;
            a = (double)peaks.get(i).x;
            
            // select b and c coefficients
            if ((modeIndex > 0) && (modeIndex < nClasses-1)){ // entry in between
                b = freq[modeIndex-1];
                c = freq[modeIndex+1];
            } else if (modeIndex == 0){ // first entry
                c = freq[modeIndex+1];
                b = (ring == true) ? freq[freq.length-1] : 0;
            } else if (modeIndex == freq.length-1) { // last entry
                b = freq[modeIndex-1];
                c = (ring == true) ? freq[0] : 0;
            }
            // calculate modes
            modes[i] = m0 + step*(a-b)/(2*a-b-c);
        }
        
        
        return modes;
    }
    
    /**
     * meanDelta function is used to calculate mean difference between
     * neighbouring elements in sorted array
     * @param values
     * @return 
     */
    static double meanDelta(double[] values){
        double dnorm = 0;
        int n = values.length-1;
        for (int i = 1; i < values.length; i++) {
            dnorm += (values[i]-values[i-1])/n; // normalized delta values
        }
        return dnorm;
    }
    
    public static int[] Histogram(double[] values, double step, boolean sort){
        if(sort){
           Arrays.sort(values); 
        }
        
        double min = values[0];
        double max = values[values.length-1];
        
        // calculate the number of classes
        int nClasses = (int) ((max-min) / step) + (int) Math.ceil((max-min) % step);
        
        int freq[] = new int[nClasses]; // frequencies array
        
        // reset all frequencies to zero
        for (int i = 0; i < freq.length; i++){
            freq[i] = 0;
        }
        
        // ceil is upper boundary for current class 
        double ceil = min + step;
                    
        int j = 0;

        // calculate frequencies
        for (int i = 0; i < freq.length; i++){
            if (j<values.length){
                while (values[j] < ceil){
                    freq[i]++;
                    j++;
                    if (j == values.length) break;
                }
            }
            ceil += step;
        }
        
        return freq;
    }
    
    public static double Min(double[] values, boolean sort){
        if (sort){
            Arrays.sort(values);
        }
        return values[0];
    }
    
    public static double Max(double[] values, boolean sort){
        if (sort){
            Arrays.sort(values);
        }
        return values[values.length-1];
    }
    
    public static double Range(double[] values, boolean sort){
        if (sort){
            Arrays.sort(values);
        }
        return Max(values, false) - Min(values, false);
    }
    
    public static double Median(double[] values, boolean sort){
        if (sort){
            Arrays.sort(values);
        }
        return values[values.length/2-1];
    }
    
    public static double Mean(double[] values){
        double mean = 0;
        for (int i = 0; i < values.length; i++) {
            mean += values[i];
        }
        
        mean /= values.length;
        return mean;
    }
    
    
    public static double MRSQ(double[] values){
        double mean = Mean(values);
        double d = 0;
        for (int i = 0; i < values.length; i++) {
            d += Math.pow(values[i]-mean, 2);
        }
        
        d /= values.length;
                
        return Math.sqrt(d);
    }
    
    public static double MRSQ(double[] values, double mean){
        double d = 0;
        for (int i = 0; i < values.length; i++) {
            d += Math.pow(values[i]-mean, 2);
        }
        
        d /= values.length;
         
        return Math.sqrt(d);
    }
    
    public static double Skewness(double[] values){
        double mean = Mean(values);
        double s3 = Math.pow(MRSQ(values, mean), 3);
        
        double m3 = 0;
        for (int i = 0; i < values.length; i++) {
            m3 += Math.pow(values[i]-mean, 3)/values.length;
        }
        
        return m3/s3;
    }
    
    public static double Kurtosis(double[] values){
        double mean = Mean(values);
        double s4 = Math.pow(MRSQ(values, mean), 4);
        
        double m4 = 0;
        for (int i = 0; i < values.length; i++) {
            m4 += Math.pow(values[i]-mean, 4)/values.length;
        }
        
        return m4/s4 - 3;
    }
    
    public static double[] ClassifyEqual(double[] values, int quantity){
        throw null;
    }
    
    public static double[] ClassifyDefined(double[] values, double step, double basevalue){
        throw null;
    }
    
    public static double[] ClassifyNatural(double[] values, int nClasses){
        //int numclass;
        int numdata = values.length;


        double[][] mat1 = new double[numdata + 1][nClasses + 1];
        double[][] mat2 = new double[numdata + 1][nClasses + 1];
        double[] st = new double[numdata];


        for (int i = 1; i <= nClasses; i++) {
            mat1[1][i] = 1;
            mat2[1][i] = 0;
            for (int j = 2; j <= numdata; j++){
                mat2[j][i] = Double.MAX_VALUE;
            }
        }
        double v = 0;
        for (int l = 2; l <= numdata; l++) {
            double s1 = 0;
            double s2 = 0;
            double w = 0;
            for (int m = 1; m <= l; m++) {
                int i3 = l - m + 1;

                double val = values[i3-1];

                s2 += val * val;
                s1 += val;

                w++;
                v = s2 - (s1 * s1) / w;
                int i4 = i3 - 1;
                if (i4 != 0) {
                    for (int j = 2; j <= nClasses; j++) {
                        if (mat2[l][j] >= (v + mat2[i4][j - 1])) {
                            mat1[l][j] = i3;
                            mat2[l][j] = v + mat2[i4][j - 1];
                        }
                    }
                }
            }
            mat1[l][1] = 1;
            mat2[l][1] = v;
        }
        int k = numdata;


        // find break items
        int[] kclass = new int[nClasses];


        kclass[nClasses - 1] = values.length - 1;


        for (int j = nClasses; j >= 2; j--) {
            System.out.println("rank = " + mat1[k][j]);
            int id =  (int) (mat1[k][j]) - 2;
            System.out.println("val = " + values[id]);

            kclass[j - 2] = id;

            k = (int) mat1[k][j] - 1;
        }
        
        double[] breaks = new double[nClasses];
        
        // find break values
        for (int i = 0; i < kclass.length; i++){
            breaks[i] = values[kclass[i]];
        }
        
        return breaks;
    }
    
    public static double[] ClassifyGeometric(double[] values, int quantity){
        throw null;
    }
    
    public static double[] ClassifyQuantile(double[] values, int quantity){
        throw null;
    }
    
    /**
     * NESTED means classification
     * @param values
     * @param quantity can be 2^x: 2, 4, 8, 16 etc.
     * @return 
     */
    public static double[] ClassifyNested(double[] values, int quantity){
        throw null;
    }
    
    /**
     * Head/Tail classification. 
     * It is a particular case of nested means classification
     * @param values
     * @param quantity
     * @return 
     */
    public static double[] ClassifyHeadtail(double[] values, int quantity){
        throw null;
    }
    
    /**
     * GetBreaks returns break values that correspond to array indexes 
     * in <code>breakitems</code> array
     * @param values
     * @param breakitems
     * @param sort
     * @return 
     */
    public static double[] getBreaks(double[] values, int breakitems[], boolean sort){
        
        // sort array if needed
        if (sort){
            Arrays.sort(values);
        }
        
        // read breaks
        double[] breaks = new double[breakitems.length];
        for (int i = 0; i < breaks.length; i++){
            breaks[i] = values[breakitems[i]];
        }
        
        return breaks;
    }
    
    /**
     * GetBreakItems returns break item indexes that correspond to break values 
     * in <code>breaks</code> array. An interesting option of this tool is that
     * it returns only those items that fall within brakes range.
     * 
     * @param values
     * @param breaks
     * @param sort
     * @return 
     */
    public static int[] getBreakItems(double[] values, double breaks[], boolean sort){
        
        // sort array if needed
        if (sort){
            Arrays.sort(values);
        }
        
        // read breaks
        ArrayList<Integer> breakItemsList = new ArrayList<>();
        
        double curBreak = breaks[0];
        
        int j = 1; // the counter of breaks
        
        for (int i = 0; i < values.length; i++){
            if(values[i] <= curBreak){
                continue;
            } else {
                breakItemsList.add(i-1);
                curBreak = breaks[j++];
            }
        }
        
        int[] breakItemsArray = new int[breakItemsList.size()];
        
        for (int i = 0; i < breakItemsArray.length; i++) {
            breakItemsArray[i] = breakItemsList.get(i);
        }
        
        return breakItemsArray;
    }
    
    /**
     * Extract range function extracts all values from <code>values</code> array 
     * that fall between <code>from</code> and <code>to</code> values
     * @param values
     * @param min
     * @param max
     * @param sort
     * @return 
     */
    public static double[] extractRange(double values[], double min, double max, boolean sort){
        if (sort){
            Arrays.sort(values);
        }
        ArrayList<Double> resultList = new ArrayList<>();
        
        // add low bound
        resultList.add(min);
        
        int i = 0;
        while(values[i] < min){
           i++;
           if (i == values.length)
               return null;
        }
        
        while(values[i] < max){
            resultList.add(values[i]);
            i++;
            if (i == values.length)
                break;
        }
        
        // add high bound
        resultList.add(max);
        
        double[] result = new double[resultList.size()];
        for (int j = 0; j < result.length; j++) {
            result[j] = resultList.get(j);
        }
        
        return result;
        
    }
}
