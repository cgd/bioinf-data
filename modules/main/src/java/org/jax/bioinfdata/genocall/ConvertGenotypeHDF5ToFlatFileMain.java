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
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jax.util.io.CommonFlatFileFormat;
import org.jax.util.io.FlatFileWriter;
import org.jax.util.io.IllegalFormatException;

import ch.systemsx.cisd.hdf5.HDF5FactoryProvider;
import ch.systemsx.cisd.hdf5.IHDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ConvertGenotypeHDF5ToFlatFileMain
{
    /**
     * the main entry point
     * @param args  command line args
     * @throws IOException
     * @throws IllegalFormatException
     */
    public static void main(String[] args) throws IllegalFormatException, IOException
    {
        // Deal with the options.
        CommandLineParser parser = new GnuParser();
        Options options = new Options();
        CommandLine commandLine = null;
        
        final Option helpOption;
        {
            helpOption = new Option("help", "Print this help and exit");
            helpOption.setRequired(false);
            options.addOption(helpOption);
        }
        
        final Option genoFileOption;
        {
            genoFileOption = new Option("genooutfile", "the genotype output flat file");
            genoFileOption.setRequired(true);
            genoFileOption.setArgs(1);
            genoFileOption.setArgName("file name");
            options.addOption(genoFileOption);
        }
        
        final Option genoOutFormatOption;
        {
            genoOutFormatOption = new Option(
                    "genooutformat",
                    "[optional] the format of the genotype file (must be \"csv\" or \"tab\")");
            genoOutFormatOption.setRequired(false);
            genoOutFormatOption.setArgs(1);
            genoOutFormatOption.setArgName("csv or tab");
            options.addOption(genoOutFormatOption);
        }
        
        final Option hdf5InputFileOption;
        {
            hdf5InputFileOption = new Option(
                    "hdf5in",
                    "the file to read HDF5 input from");
            hdf5InputFileOption.setRequired(true);
            hdf5InputFileOption.setArgs(1);
            hdf5InputFileOption.setArgName("file name");
            options.addOption(hdf5InputFileOption);
        }
        
        try
        {
            commandLine = parser.parse(options, args);
            
            // See if we just need to print the help options.
            if(commandLine.hasOption(helpOption.getOpt()))
            {
                HelpFormatter helpFormatter = new HelpFormatter();
                helpFormatter.printHelp("hdf5toff", options);
            }
            else
            {
                final String genoFileName = commandLine.getOptionValue(genoFileOption.getOpt());
                final String genoOutFmtStr = commandLine.getOptionValue(genoOutFormatOption.getOpt());
                final String hdf5InFileName = commandLine.getOptionValue(hdf5InputFileOption.getOpt());
                
                IHDF5Factory hdf5Fac = HDF5FactoryProvider.get();
                IHDF5Reader hdf5Reader = hdf5Fac.openForReading(new File(hdf5InFileName));
                HDF5GenotypeCallMatrix hdf5GenoMatrix = new HDF5GenotypeCallMatrix(hdf5Reader);
                
                final FlatFileWriter genoFFW;
                if(genoOutFmtStr == null || genoOutFmtStr.trim().toLowerCase().equals("csv"))
                {
                    genoFFW = new FlatFileWriter(
                            new FileWriter(genoFileName),
                            CommonFlatFileFormat.CSV_UNIX);
                }
                else if(genoOutFmtStr.trim().toLowerCase().equals("tab"))
                {
                    genoFFW = new FlatFileWriter(
                            new FileWriter(genoFileName),
                            CommonFlatFileFormat.TAB_DELIMITED_UNIX);
                }
                else
                {
                    throw new ParseException("geno file input format must be \"tab\" or \"csv\"");
                }
                
                GenotypesFlatFile gff = new GenotypesFlatFile();
                gff.writeGenoCallMatrix(hdf5GenoMatrix, genoFFW);
                hdf5Reader.close();
                genoFFW.close();
            }
        }
        catch(ParseException ex)
        {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp("hdf5toff", options);
            
            System.exit(-1);
        }
    }
}
