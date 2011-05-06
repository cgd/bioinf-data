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
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jax.util.io.CommonFlatFileFormat;
import org.jax.util.io.FlatFileReader;
import org.jax.util.io.IllegalFormatException;

import ch.systemsx.cisd.hdf5.HDF5FactoryProvider;
import ch.systemsx.cisd.hdf5.IHDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ConvertAlchemyCallsToHDF5Main
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
            genoFileOption = new Option("alchemygenos", "the genotype calls output from alchemy");
            genoFileOption.setRequired(true);
            genoFileOption.setArgs(1);
            genoFileOption.setArgName("file name");
            options.addOption(genoFileOption);
        }
        
        final Option bpBuildIDOption;
        {
            bpBuildIDOption = new Option(
                    "bpbuild",
                    "[optional] BP position build identifier (Eg: \"NCBI Build 37\")");
            bpBuildIDOption.setRequired(false);
            bpBuildIDOption.setArgs(1);
            bpBuildIDOption.setArgName("build identifier");
            options.addOption(bpBuildIDOption);
        }
        
        final Option outputFileOption;
        {
            outputFileOption = new Option(
                    "hdf5out",
                    "the file to write HDF5 output to");
            outputFileOption.setRequired(true);
            outputFileOption.setArgs(1);
            outputFileOption.setArgName("file name");
            options.addOption(outputFileOption);
        }
        
        try
        {
            commandLine = parser.parse(options, args);
            
            // See if we just need to print the help options.
            if(commandLine.hasOption(helpOption.getOpt()))
            {
                HelpFormatter helpFormatter = new HelpFormatter();
                helpFormatter.printHelp("alchtohdf5", options);
            }
            else
            {
                final String genoFileName = commandLine.getOptionValue(genoFileOption.getOpt());
                final String bpBuildIDStr = commandLine.getOptionValue(bpBuildIDOption.getOpt());
                final String outFileName = commandLine.getOptionValue(outputFileOption.getOpt());
                
                final FlatFileReader genoFFR = new FlatFileReader(
                        new FileReader(genoFileName),
                        CommonFlatFileFormat.TAB_DELIMITED_UNIX);
                
                GenotypesFlatFile gff = new GenotypesFlatFile();
                GenotypeCallMatrix genoMat = gff.readAlchemyGenoCalls(genoFFR, bpBuildIDStr);
                genoFFR.close();
                
                GenotypesHDF5 ghdf5 = new GenotypesHDF5();
                IHDF5Factory hdf5Fac = HDF5FactoryProvider.get();
                File hdf5File = new File(outFileName);
                if(hdf5File.exists())
                {
                    if(!hdf5File.delete())
                    {
                        throw new IOException(
                                "failed to overwrite \"" + outFileName + "\"");
                    }
                }
                IHDF5Writer hdf5Writer = hdf5Fac.open(hdf5File);
                ghdf5.writeGenoCallMatrix(genoMat, hdf5Writer);
                hdf5Writer.close();
            }
        }
        catch(ParseException ex)
        {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp("fftohdf5", options);
            
            System.exit(-1);
        }
    }
}
