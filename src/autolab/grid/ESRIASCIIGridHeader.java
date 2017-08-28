/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package autolab.grid;

import java.io.IOException;
import java.util.Scanner;

/**
 *
 * @author tsamsonov
 */
public class ESRIASCIIGridHeader {

        public int cols = 0;
        public int rows = 0;
        public double west = Double.NaN;
        public double south = Double.NaN;
        public double cellSize = Double.NaN;
        public float noDataValue = Float.NaN;

        /*
         * returns whether valid values have been found
         */
        protected boolean isValid() {

            return (cols > 0 
                    && rows > 0
                    && cellSize > 0
                    && !Double.isNaN(west)
                    && !Double.isNaN(south));
            // noDataValue is optional
        }

        /**
         * Reads cols, rows, west, south, cellSize and noDataValue from the header.
         * Throws an exception if this is not a valid header.
         * @param scanner Scanner must be initialized to use dot as decimal separator.
         * @throws IOException
         */
        public void readHeader(Scanner scanner) throws IOException {

            cols = rows = 0;
            west = south = cellSize = Double.NaN;
            noDataValue = Float.NaN;

            while (scanner.hasNext()) {

                if (scanner.hasNextDouble()) {
                    // next line starts with number, must be grid
                    break;
                }

                String str = scanner.next().trim().toLowerCase();
                if (str.equals("ncols")) {
                    this.cols = scanner.nextInt();
                } else if (str.equals("nrows")) {
                    this.rows = scanner.nextInt();
                } else if (str.equals("xllcenter") || str.equals("xllcorner")) {
                    this.west = scanner.nextDouble();
                } else if (str.equals("yllcenter") || str.equals("yllcorner")) {
                    this.south = scanner.nextDouble();
                } else if (str.equals("cellsize")) {
                    this.cellSize = scanner.nextDouble();
                } else if (str.startsWith("nodata")) {
                    this.noDataValue = scanner.nextFloat();
                } else {

                    // make sure the line starts with a number
                    if (!scanner.hasNextDouble()) {
                        throw new IOException();
                    }

                    // done reading the header
                    break;
                }
            }

            if (!isValid()) {
                throw new IOException();
            }
        }
    }

