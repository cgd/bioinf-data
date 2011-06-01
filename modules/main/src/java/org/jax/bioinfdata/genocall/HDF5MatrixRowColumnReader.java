/*
 * Copyright (c) 2010 The Jackson Laboratory
 * 
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.bioinfdata.genocall;

import java.io.File;

import org.jax.util.datastructure.SequenceUtilities;

import ch.systemsx.cisd.hdf5.HDF5FactoryProvider;
import ch.systemsx.cisd.hdf5.IHDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

/**
 * A class that provides some convenience methods for getting matrix rows and
 * columns (rows are considered dimension 0 and columns are dimension 1)
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HDF5MatrixRowColumnReader
{
    private final IHDF5Reader reader;
    private final String objectPath;
    private final long rowCount;
    private final long columnCount;
    
    /**
     * Constructor
     * @param reader
     *          the HDF5 reader
     * @param objectPath
     *          the path to the matrix object
     */
    public HDF5MatrixRowColumnReader(IHDF5Reader reader, String objectPath)
    {
        super();
        this.reader = reader;
        this.objectPath = objectPath;
        long[] dimensions = reader.getDataSetInformation(objectPath).getDimensions();
        if(dimensions.length == 2)
        {
            this.rowCount = dimensions[0];
            this.columnCount = dimensions[1];
        }
        else
        {
            throw new IllegalArgumentException(
                    "The given object was expected to have two dimensions but " +
                    "instead it has " + dimensions.length);
        }
    }
    
    /**
     * Getter for the row count
     * @return the rowCount
     */
    public long getRowCount()
    {
        return this.rowCount;
    }
    
    /**
     * Getter for the column count
     * @return the column count
     */
    public long getColumnCount()
    {
        return this.columnCount;
    }
    
    /**
     * Get bytes for the given row
     * @param index
     *          the row index
     * @return
     *          the row bytes
     */
    public byte[] getRowBytes(long index)
    {
        byte[][] block = this.reader.readByteMatrixBlock(
                this.objectPath,
                1, (int)this.columnCount,
                index, 0);
        return block[0];
    }
    
    /**
     * Get bytes for the column
     * @param index
     *          the column index
     * @return
     *          the bytes
     */
    public byte[] getColumnBytes(long index)
    {
        byte[][] block = this.reader.readByteMatrixBlock(
                this.objectPath,
                (int)this.rowCount, 1,
                0, index);
        byte[] col = new byte[(int)this.rowCount];
        for(int row = 0; row < this.rowCount; row++)
        {
            col[row] = block[row][0];
        }
        return col;
    }
    
    /**
     * Just for test purposes
     * @param args
     */
    public static void main(String[] args)
    {
        IHDF5Factory hdf5Fac = HDF5FactoryProvider.get();
        IHDF5Reader hdf5Reader = hdf5Fac.openForReading(new File(args[0]));
        
        HDF5MatrixRowColumnReader rowColReader = new HDF5MatrixRowColumnReader(
                hdf5Reader,
                args[1]);
        System.out.println("printing by row");
        for(int rowIndex = 0; rowIndex < rowColReader.getRowCount(); rowIndex++)
        {
            byte[] row = rowColReader.getRowBytes(rowIndex);
            System.out.println(SequenceUtilities.toString(
                    SequenceUtilities.toByteList(row), ","));
        }
        System.out.println("printing by col");
        for(int colIndex = 0; colIndex < rowColReader.getColumnCount(); colIndex++)
        {
            byte[] col = rowColReader.getColumnBytes(colIndex);
            System.out.println(SequenceUtilities.toString(
                    SequenceUtilities.toByteList(col), ","));
        }
    }
}
