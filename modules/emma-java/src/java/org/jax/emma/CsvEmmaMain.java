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

package org.jax.emma;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jax.bioinfdata.genocall.GenotypeCallMatrix;
import org.jax.bioinfdata.genocall.GenotypesFlatFile;
import org.jax.util.io.CommonFlatFileFormat;
import org.jax.util.io.FlatFileReader;
import org.jax.util.io.IllegalFormatException;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CsvEmmaMain
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
            genoFileOption = new Option("genofile", "the genotype CSV file");
            genoFileOption.setRequired(true);
            genoFileOption.setArgs(1);
            genoFileOption.setArgName("file name");
            options.addOption(genoFileOption);
        }
        
        final Option aAlleleOption;
        {
            aAlleleOption = new Option("aallelecol", "the A allele column # (one-based index)");
            aAlleleOption.setRequired(true);
            aAlleleOption.setArgs(1);
            aAlleleOption.setArgName("column number");
            options.addOption(aAlleleOption);
        }
        
        final Option bAlleleOption;
        {
            bAlleleOption = new Option("ballelecol", "the B allele column # (one-based index)");
            bAlleleOption.setRequired(true);
            bAlleleOption.setArgs(1);
            bAlleleOption.setArgName("column number");
            options.addOption(bAlleleOption);
        }
        
        final Option firstGenoColumnOption;
        {
            firstGenoColumnOption = new Option("firstgenocol", "the first genotype column # (one-based index)");
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
        
        final Option phenoFileOption;
        {
            phenoFileOption = new Option("phenofile", "the phenotype file");
            phenoFileOption.setRequired(true);
            phenoFileOption.setArgs(1);
            phenoFileOption.setArgName("file name");
            options.addOption(phenoFileOption);
        }
        
        final Option phenoNameOption;
        {
            phenoNameOption = new Option("phenoname", "[optional] the name of the phenotype to scan");
            phenoNameOption.setRequired(false);
            phenoNameOption.setArgs(1);
            phenoNameOption.setArgName("name");
            options.addOption(phenoNameOption);
        }
        
        final Option sexOption;
        {
            sexOption = new Option(
                    "sex",
                    "[optional] filter phenotypes by sex. " +
                    "agnostic is the default value.");
            sexOption.setRequired(false);
            sexOption.setArgs(1);
            sexOption.setArgName("agnostic/female/male");
            options.addOption(sexOption);
        }
        
        final Option outputFileOption;
        {
            outputFileOption = new Option("out", "the file to write scan output to");
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
                helpFormatter.printHelp("emmascan", options);
            }
            else
            {
                final String genoFileName = commandLine.getOptionValue(genoFileOption.getOpt());
                final String aColStr = commandLine.getOptionValue(aAlleleOption.getOpt());
                final String bColStr = commandLine.getOptionValue(bAlleleOption.getOpt());
                final String fstGenoColStr = commandLine.getOptionValue(firstGenoColumnOption.getOpt());
                final String lstGenoColStr = commandLine.getOptionValue(lastGenoColumnOption.getOpt());
                final String phenoFileName = commandLine.getOptionValue(phenoFileOption.getOpt());
                final String phenotype = commandLine.getOptionValue(phenoNameOption.getOpt());
                final String sexStr = commandLine.getOptionValue(sexOption.getOpt());
                final String outFileName = commandLine.getOptionValue(outputFileOption.getOpt());
                
                final SexFilter sexToScan;
                if(sexStr == null || sexStr.toLowerCase().equals("agnostic"))
                {
                    sexToScan = SexFilter.AGNOSTIC;
                }
                else if(sexStr.toLowerCase().equals("female"))
                {
                    sexToScan = SexFilter.ALLOW_FEMALE;
                }
                else if(sexStr.toLowerCase().equals("male"))
                {
                    sexToScan = SexFilter.ALLOW_MALE;
                }
                else
                {
                    throw new ParseException("sex option cannot be: " + sexStr);
                }
                
                FlatFileReader ffr = new FlatFileReader(
                        new BufferedReader(new FileReader(genoFileName)),
                        CommonFlatFileFormat.CSV_UNIX);
                GenotypeCallMatrix genoMat = GenotypesFlatFile.readGenoCallMatrix(
                        new FlatFileReader[] {ffr},
                        Integer.parseInt(aColStr.trim()) - 1,
                        Integer.parseInt(bColStr.trim()) - 1,
                        -1,
                        -1,
                        -1,
                        null,
                        Integer.parseInt(fstGenoColStr.trim()) - 1,
                        lstGenoColStr == null ? -1 : Integer.parseInt(lstGenoColStr.trim()));
                ffr.close();
                
                EMMAAssociationTest emmaTest = new EMMAAssociationTest();
                double[] scanResults = emmaTest.emmaScan(
                        genoMat,
                        phenoFileName,
                        phenotype,
                        sexToScan);
                
                PrintStream out = new PrintStream(outFileName);
                out.println("pValue");
                for(int i = 0; i < scanResults.length; i++)
                {
                    out.println(scanResults[i]);
                }
                out.close();
            }
        }
        catch(ParseException ex)
        {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp("emmascan", options);
            
            System.exit(-1);
        }
    }
}
