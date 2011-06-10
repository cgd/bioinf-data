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
import org.jax.util.io.FlatFileFormat;
import org.jax.util.io.FlatFileReader;
import org.jax.util.io.IllegalFormatException;

import ch.systemsx.cisd.hdf5.HDF5FactoryProvider;
import ch.systemsx.cisd.hdf5.IHDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ConvertGenotypeFlatFileToHDF5Main
{
    private static final int argToIntMinus1(String arg)
    {
        return arg == null ? -1 : Integer.parseInt(arg.trim()) - 1;
    }

    private static final int argToInt(String arg)
    {
        return arg == null ? -1 : Integer.parseInt(arg.trim());
    }

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
            genoFileOption = new Option("genoinfiles", "the genotype input flat file");
            genoFileOption.setRequired(true);
            genoFileOption.setArgs(Integer.MAX_VALUE);
            genoFileOption.setArgName("file name");
            options.addOption(genoFileOption);
        }
        
        final Option genoInFormatOption;
        {
            genoInFormatOption = new Option(
                    "genoinformat",
                    "[optional] the format of the genotype file (must be \"csv\"" +
                    " or \"tab\". \"csv\" is the default)");
            genoInFormatOption.setRequired(false);
            genoInFormatOption.setArgs(1);
            genoInFormatOption.setArgName("csv or tab");
            options.addOption(genoInFormatOption);
        }
        
        final Option aAlleleOption;
        {
            aAlleleOption = new Option(
                    "aallelecol",
                    "[optional] the A allele column # (one-based index). If " +
                    "no A allele column is given then allele codes must be " +
                    "used in place of nucleotide values where 1 = A allele, " +
                    "2 = B allele, 3 = Heterozygous, -1 = No Call");
            aAlleleOption.setRequired(false);
            aAlleleOption.setArgs(1);
            aAlleleOption.setArgName("column #");
            options.addOption(aAlleleOption);
        }
        
        final Option bAlleleOption;
        {
            bAlleleOption = new Option(
                    "ballelecol",
                    "[optional] the B allele column # (one-based index). If " +
                    "no B allele column is given then allele codes must be " +
                    "used in place of nucleotide values where 1 = A allele, " +
                    "2 = B allele, 3 = Heterozygous, -1 = No Call");
            bAlleleOption.setRequired(false);
            bAlleleOption.setArgs(1);
            bAlleleOption.setArgName("column #");
            options.addOption(bAlleleOption);
        }
        
        final Option snpIdColumnOption;
        {
            snpIdColumnOption = new Option(
                    "snpcol",
                    "[optional] the SNP ID column # (one-based index)");
            snpIdColumnOption.setRequired(false);
            snpIdColumnOption.setArgs(1);
            snpIdColumnOption.setArgName("column #");
            options.addOption(snpIdColumnOption);
        }
        
        final Option chromosomeColumnOption;
        {
            chromosomeColumnOption = new Option(
                    "chrcol",
                    "[optional] the chromosome ID column # (one-based index)");
            chromosomeColumnOption.setRequired(false);
            chromosomeColumnOption.setArgs(1);
            chromosomeColumnOption.setArgName("column #");
            options.addOption(chromosomeColumnOption);
        }
        
        final Option bpPositionColumnOption;
        {
            bpPositionColumnOption = new Option(
                    "poscol",
                    "[optional] base-pair position column # (one-based index)");
            bpPositionColumnOption.setRequired(false);
            bpPositionColumnOption.setArgs(1);
            bpPositionColumnOption.setArgName("column #");
            options.addOption(bpPositionColumnOption);
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
        
        final Option firstGenoColumnOption;
        {
            firstGenoColumnOption = new Option(
                    "firstgenocol",
                    "the first genotype column # (one-based index)");
            firstGenoColumnOption.setRequired(true);
            firstGenoColumnOption.setArgs(1);
            firstGenoColumnOption.setArgName("column #");
            options.addOption(firstGenoColumnOption);
        }
        
        final Option lastGenoColumnOption;
        {
            lastGenoColumnOption = new Option(
                    "lastgenocol",
                    "[optional] the last genotype column # (one-based index). " +
                    "The default behavior is to assume that all columns after " +
                    "the first genotype column are genotype columns.");
            lastGenoColumnOption.setRequired(false);
            lastGenoColumnOption.setArgs(1);
            lastGenoColumnOption.setArgName("column #");
            options.addOption(lastGenoColumnOption);
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
                helpFormatter.printHelp("fftohdf5", options);
            }
            else
            {
                final String[] genoFileNames = commandLine.getOptionValues(genoFileOption.getOpt());
                final String genoInFmtStr = commandLine.getOptionValue(genoInFormatOption.getOpt());
                final String aColStr = commandLine.getOptionValue(aAlleleOption.getOpt());
                final String bColStr = commandLine.getOptionValue(bAlleleOption.getOpt());
                final String snpColStr = commandLine.getOptionValue(snpIdColumnOption.getOpt());
                final String chrColStr = commandLine.getOptionValue(chromosomeColumnOption.getOpt());
                final String bpPosStr = commandLine.getOptionValue(bpPositionColumnOption.getOpt());
                final String bpBuildIDStr = commandLine.getOptionValue(bpBuildIDOption.getOpt());
                final String fstGenoColStr = commandLine.getOptionValue(firstGenoColumnOption.getOpt());
                final String lstGenoColStr = commandLine.getOptionValue(lastGenoColumnOption.getOpt());
                final String outFileName = commandLine.getOptionValue(outputFileOption.getOpt());
                
                final FlatFileFormat ffFormat;
                if(genoInFmtStr == null || genoInFmtStr.trim().toLowerCase().equals("csv"))
                {
                    ffFormat = CommonFlatFileFormat.CSV_UNIX;
                }
                else if(genoInFmtStr.trim().toLowerCase().equals("tab"))
                {
                    ffFormat = CommonFlatFileFormat.TAB_DELIMITED_UNIX;
                }
                else
                {
                    throw new ParseException("geno file input format must be \"tab\" or \"csv\"");
                }
                
                final FlatFileReader[] genoFFRs = new FlatFileReader[genoFileNames.length];
                for(int i = 0; i < genoFFRs.length; i++)
                {
                    genoFFRs[i] = new FlatFileReader(
                            new FileReader(genoFileNames[i]),
                            ffFormat);
                }
                
                GenotypesFlatFile gff = new GenotypesFlatFile();
                GenotypeCallMatrix genoMat = gff.readGenoCallMatrix(
                        genoFFRs,
                        argToIntMinus1(aColStr),
                        argToIntMinus1(bColStr),
                        argToIntMinus1(snpColStr),
                        argToIntMinus1(chrColStr),
                        argToIntMinus1(bpPosStr),
                        bpBuildIDStr,
                        argToIntMinus1(fstGenoColStr),
                        argToInt(lstGenoColStr));
                
                for(int i = 0; i < genoFFRs.length; i++)
                {
                    genoFFRs[i].close();
                }
                
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
                HDF5GenotypeCallMatrix hdf5GenoMat = new HDF5GenotypeCallMatrix(hdf5Writer);
                AbstractGenotypeCallMatrix.copyGenoMatrix(genoMat, hdf5GenoMat);
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
