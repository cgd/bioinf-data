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

package org.jax.bioinfdata.phylogeny;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jax.bioinfdata.genocall.AbstractGenotypeCallMatrix;
import org.jax.bioinfdata.genocall.HDF5GenotypeCallMatrix;
import org.jax.util.io.CommonFlatFileFormat;
import org.jax.util.io.FlatFileWriter;

import ch.systemsx.cisd.hdf5.HDF5FactoryProvider;
import ch.systemsx.cisd.hdf5.IHDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class MaxKPhylogenyMain
{
    /**
     * The main entry point
     * @param args function arguments
     * @throws IOException if we have a problem reading/writing data
     * @throws NoValidPhylogenyException
     */
    public static void main(String[] args) throws IOException, NoValidPhylogenyException
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
        
        final Option hdf5InputFileOption;
        {
            hdf5InputFileOption = new Option(
                    "hdf5in",
                    "the HDF5 input file containing the genotype matrix");
            hdf5InputFileOption.setRequired(true);
            hdf5InputFileOption.setArgs(1);
            hdf5InputFileOption.setArgName("file name");
            options.addOption(hdf5InputFileOption);
        }
        
        final Option csvOutputFileOption;
        {
            csvOutputFileOption = new Option(
                    "csvout",
                    "the CSV output file with max-k intervals and their " +
                    "corresponding perfect phylogenies in newick format");
            csvOutputFileOption.setRequired(true);
            csvOutputFileOption.setArgs(1);
            csvOutputFileOption.setArgName("file name");
            options.addOption(csvOutputFileOption);
        }
        
        try
        {
            commandLine = parser.parse(options, args);
            
            // See if we just need to print the help options.
            if(commandLine.hasOption(helpOption.getOpt()))
            {
                HelpFormatter helpFormatter = new HelpFormatter();
                helpFormatter.printHelp("maxkphylo", options);
            }
            else
            {
                final String hdf5InFileName = commandLine.getOptionValue(hdf5InputFileOption.getOpt());
                final String csvOutFileName = commandLine.getOptionValue(csvOutputFileOption.getOpt());
                
                IHDF5Factory hdf5Fac = HDF5FactoryProvider.get();
                IHDF5Reader hdf5Reader = hdf5Fac.openForReading(new File(hdf5InFileName));
                HDF5GenotypeCallMatrix hdf5GenoMatrix = new HDF5GenotypeCallMatrix(hdf5Reader);
                
                // write the header row
                FlatFileWriter ffw = new FlatFileWriter(
                        //new OutputStreamWriter(System.out),
                        new FileWriter(csvOutFileName),
                        CommonFlatFileFormat.CSV_UNIX);
                ffw.writeRow(new String[] {
                        "chrID",
                        "bpStartPosition",
                        "bpEndPosition",
                        "newickPerfectPhylogeny"});
                
                // scan the chromosomes in order
                for(AbstractGenotypeCallMatrix currChrView : hdf5GenoMatrix.getChromosomeViews())
                {
                    List<IndexedSnpInterval> maxKScanResult = IntervalScanner.maxKScan(currChrView);
                    List<PhylogenyTreeNode> phyloScanResult = PhylogenyScanner.inferPerfectPhylogenies(
                            currChrView,
                            maxKScanResult);
                    assert maxKScanResult.size() == phyloScanResult.size();
                    
                    String[] chrIDArray = currChrView.getChrIDs();
                    long[] posArray = currChrView.getBpPositions();
                    for(int i = 0; i < maxKScanResult.size(); i++)
                    {
                        IndexedSnpInterval isi = maxKScanResult.get(i);
                        PhylogenyTreeNode phylo = phyloScanResult.get(i);
                        ffw.writeRow(new String[] {
                                chrIDArray[isi.getStartIndex()],
                                Long.toString(posArray[isi.getStartIndex()]),
                                Long.toString(posArray[isi.getEndIndex()]),
                                phylo.toNewickFormat()});
                    }
                }
                ffw.flush();
                ffw.close();
            }
        }
        catch(ParseException ex)
        {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp("maxkphylo", options);
            
            System.exit(-1);
        }
    }
}
