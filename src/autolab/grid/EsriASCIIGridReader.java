/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.grid;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Locale;
import java.util.Scanner;
import java.util.StringTokenizer;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridCoverageFactory;
import org.geotools.geometry.jts.ReferencedEnvelope;

/**
 * A simple reader for Esri ASCII grid files.
 * Adapted from Bernhard Jenny's (jenny@oregonstate.edu) code 
 * for GetTools Coverage2D format instead of GeoGrid.
 * @author jenny, tsamsonov
 */
public class EsriASCIIGridReader {
    
    /**
     * Returns whether a scanner references valid data that can be read.
     * @param scanner
     * @return
     * @throws IOException
     */
    public static boolean canRead(Scanner scanner) {
        try {
            ESRIASCIIGridHeader header = new ESRIASCIIGridHeader();
            header.readHeader(scanner);
            return header.isValid();
        } catch (Exception exc) {
            return false;
        }
    }

    public static boolean canRead(String filePath) {
        Scanner scanner = null;
        try {
            scanner = createUSScanner(new FileInputStream(filePath));
            return EsriASCIIGridReader.canRead(scanner);
        } catch (Exception exc) {
            return false;
        } finally {
            if (scanner != null) {
                try {
                    scanner.close();
                } catch (Throwable exc) {
                }
            }
        }
    }

    /** Read a Grid from a file in ESRI ASCII format.
     * @param fileName The path to the file to be read.
     * @param progress A WorkerProgress to inform about the progress.
     * @return The read grid.
     */
    public static GridCoverage2D read(File file) throws java.io.IOException {
        FileInputStream fis = new FileInputStream(file.getAbsolutePath());
        GridCoverage2D grid = EsriASCIIGridReader.read(fis);
        String name = file.getName();
        return grid;
    }
    
    /** Read a Grid from a file in ESRI ASCII format to 2-dimension array
     * @param fileName The path to the file to be read.
     * @param progress A WorkerProgress to inform about the progress.
     * @return The read grid.
     */
    public static float[][] readMatrix(File file) throws java.io.IOException {
        FileInputStream fis = new FileInputStream(file.getAbsolutePath());
        float[][] grid = EsriASCIIGridReader.readMatrix(fis);
        String name = file.getName();
        return grid;
    }
    
    /** Read a Grid from a file in ESRI ASCII format to 2-dimension array
     * @param fileName The path to the file to be read.
     * @param progress A WorkerProgress to inform about the progress.
     * @return The read grid.
     */
    public static double[][] readDoubleMatrix(File file) throws java.io.IOException {
        FileInputStream fis = new FileInputStream(file.getAbsolutePath());
        double[][] grid = EsriASCIIGridReader.readDoubleMatrix(fis);
        String name = file.getName();
        return grid;
    }

    /** Read a Grid from a stream in ESRI ASCII format.
     * @param is The stream to read from. The stream is closed at the end.
     * @param progress A WorkerProgress to inform about the progress.
     * @return The read grid.
     */
    public static GridCoverage2D read(InputStream input)
            throws IOException {


        Scanner scanner = createUSScanner(new BufferedInputStream(input));       
        try {
            ESRIASCIIGridHeader header = new ESRIASCIIGridHeader();
            header.readHeader(scanner);

            // Create envelope
            ReferencedEnvelope envelope =  new ReferencedEnvelope(
                    header.west, header.west + header.cols * header.cellSize,
                    header.south, header.south + header.rows * header.cellSize,
                    null
            );
            
            float[][] values = new float[header.rows][header.cols];

            // use legacy StringTokenizer, which is considerably faster than
            // the Scanner class, which uses regular expressions.
            StringTokenizer tokenizer = new StringTokenizer(scanner.nextLine(), " ");

            // read grid values. Rows are ordered top to bottom.
            for (int row = 0; row < header.rows; row++) {
                // read one row
                for (int col = 0; col < header.cols; col++) {

                    // a logical row in the grid does not necesseraly correspond
                    // to a line in the file!
                    if (!tokenizer.hasMoreTokens()) {
                        tokenizer = new StringTokenizer(scanner.nextLine(), " ");
                    }
                    final float v = Float.parseFloat(tokenizer.nextToken());
                    if (v == header.noDataValue || Float.isNaN(v)) {
                        values[row][col] = Float.NaN;
                    } else {
                        values[row][col] = v;
                    }
                }
            }
            GridCoverageFactory gcf = new GridCoverageFactory();
            GridCoverage2D grid = gcf.create("Yaya", values, envelope);
            
            return grid;
            
        } finally {
            try {
                // this closes the input stream
                scanner.close();
            } catch (Exception exc) {
            }
        }

    }
    
    /** Read a Grid from a stream in ESRI ASCII format to float[][] array
     * @param is The stream to read from. The stream is closed at the end.
     * @param progress A WorkerProgress to inform about the progress.
     * @return The read grid.
     */
    public static float[][] readMatrix(InputStream input)
            throws IOException {


        Scanner scanner = createUSScanner(new BufferedInputStream(input));       
        try {
            ESRIASCIIGridHeader header = new ESRIASCIIGridHeader();
            header.readHeader(scanner);

            // Create envelope
            ReferencedEnvelope envelope =  new ReferencedEnvelope(
                    header.west, header.west + header.cols * header.cellSize,
                    header.south, header.south + header.rows * header.cellSize,
                    null
            );
            
            float[][] values = new float[header.rows][header.cols];

            // use legacy StringTokenizer, which is considerably faster than
            // the Scanner class, which uses regular expressions.
            StringTokenizer tokenizer = new StringTokenizer(scanner.nextLine(), " ");

            // read grid values. Rows are ordered top to bottom.
            for (int row = 0; row < header.rows; row++) {
                // read one row
                for (int col = 0; col < header.cols; col++) {

                    // a logical row in the grid does not necesseraly correspond
                    // to a line in the file!
                    if (!tokenizer.hasMoreTokens()) {
                        tokenizer = new StringTokenizer(scanner.nextLine(), " ");
                    }
                    final float v = Float.parseFloat(tokenizer.nextToken());
                    if (v == header.noDataValue || Float.isNaN(v)) {
                        values[row][col] = Float.NaN;
                    } else {
                        values[row][col] = v;
                    }
                }
            }
//            GridCoverageFactory gcf = new GridCoverageFactory();
//            GridCoverage2D grid = gcf.create("Yaya", values, envelope);
            
            return values;
            
        } finally {
            try {
                // this closes the input stream
                scanner.close();
            } catch (Exception exc) {
            }
        }

    }
    
    public static double[][] readDoubleMatrix(InputStream input)
            throws IOException {


        Scanner scanner = createUSScanner(new BufferedInputStream(input));       
        try {
            ESRIASCIIGridHeader header = new ESRIASCIIGridHeader();
            header.readHeader(scanner);

            // Create envelope
            ReferencedEnvelope envelope =  new ReferencedEnvelope(
                    header.west, header.west + header.cols * header.cellSize,
                    header.south, header.south + header.rows * header.cellSize,
                    null
            );
            
            double[][] values = new double[header.rows][header.cols];

            // use legacy StringTokenizer, which is considerably faster than
            // the Scanner class, which uses regular expressions.
            StringTokenizer tokenizer = new StringTokenizer(scanner.nextLine(), " ");

            // read grid values. Rows are ordered top to bottom.
            for (int row = 0; row < header.rows; row++) {
                // read one row
                for (int col = 0; col < header.cols; col++) {

                    // a logical row in the grid does not necesseraly correspond
                    // to a line in the file!
                    if (!tokenizer.hasMoreTokens()) {
                        tokenizer = new StringTokenizer(scanner.nextLine(), " ");
                    }
                    final double v = Double.parseDouble(tokenizer.nextToken());
                    if (v == header.noDataValue || Double.isNaN(v)) {
                        values[row][col] = Float.NaN;
                    } else {
                        values[row][col] = v;
                    }
                }
            }
//            GridCoverageFactory gcf = new GridCoverageFactory();
//            GridCoverage2D grid = gcf.create("Yaya", values, envelope);
            
            return values;
            
        } finally {
            try {
                // this closes the input stream
                scanner.close();
            } catch (Exception exc) {
            }
        }

    }

    /**
     * Creates a scanner for ASCII text with a period as decimal separator.
     * @param is
     * @return
     * @throws FileNotFoundException
     */
    public static Scanner createUSScanner(InputStream is) throws FileNotFoundException {
        Scanner scanner = new Scanner(is, "US-ASCII");
        scanner.useLocale(Locale.US);
        return scanner;
    }
}
