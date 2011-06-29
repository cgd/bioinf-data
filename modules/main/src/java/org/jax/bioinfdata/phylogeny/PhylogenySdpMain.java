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

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jax.util.io.CommonFlatFileFormat;
import org.jax.util.io.FlatFileReader;
import org.jax.util.io.FlatFileWriter;
import org.jax.util.io.IllegalFormatException;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenySdpMain
{
    /**
     * Main entry point for reading in newick trees at specified intervals and
     * turning them into SDPs
     * @param args
     * @throws IOException 
     * @throws IllegalFormatException 
     */
    public static void main(String[] args) throws IOException, IllegalFormatException
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
        
        final Option inputFileOption;
        {
            inputFileOption = new Option(
                    "phylocsvin",
                    "This is the input " +
                    "CSV file which must have a header row allong with " +
                    "the following four columns in order: chromosome ID, " +
                    "interval start in base pairs, interval end in base pairs, " +
                    "newick formatted phylogeny tree");
            inputFileOption.setRequired(true);
            inputFileOption.setArgs(1);
            inputFileOption.setArgName("file name");
            options.addOption(inputFileOption);
        }
        
        final Option minMinorStrainCountOption;
        {
            minMinorStrainCountOption = new Option(
                    "minorallelecountthreshold",
                    "this option specifies the minimum minor allele count that " +
                    "an SDP must have. All SDPs that fall below this threshold " +
                    "will be filtered from the output.");
            minMinorStrainCountOption.setRequired(true);
            minMinorStrainCountOption.setArgs(1);
            minMinorStrainCountOption.setArgName("min SDP count");
            options.addOption(minMinorStrainCountOption);
        }
        
        final Option outputFileOption;
        {
            outputFileOption = new Option(
                    "sdpcsvout",
                    "the output CSV file");
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
                helpFormatter.printHelp("phylosdp", options);
            }
            else
            {
                final String inFileName = commandLine.getOptionValue(inputFileOption.getOpt());
                final String outFileName = commandLine.getOptionValue(outputFileOption.getOpt());
                final String minMinorStrainCountStr = commandLine.getOptionValue(
                        minMinorStrainCountOption.getOpt());
                final int minMinorStrainCount = Integer.parseInt(minMinorStrainCountStr);
                FlatFileReader ffr = new FlatFileReader(
                        new FileReader(inFileName),
                        CommonFlatFileFormat.CSV_UNIX);
                String[] inHeader = ffr.readRow();
                if(inHeader == null)
                {
                    throw new IOException("the input file is empty");
                }
                else if(inHeader.length != 4)
                {
                    throw new IllegalFormatException(
                            "expected the input to have 4 columns but found " +
                            inHeader.length + " columns");
                }
                else
                {
                    Map<BitSet, List<GenomeInterval>> sdpMap =
                        new HashMap<BitSet, List<GenomeInterval>>();
                    
                    ArrayList<String> strainNames = null;
                    String[] outHeader = null;
                    String[] currInRow = null;
                    while((currInRow = ffr.readRow()) != null)
                    {
                        GenomeInterval genoInt = new GenomeInterval(
                                currInRow[0],
                                Long.parseLong(currInRow[1]),
                                Long.parseLong(currInRow[2]));
                        PhylogenyTreeNode phylo = PhylogenyTreeNode.fromNewickFormat(currInRow[3]);
                        if(strainNames == null)
                        {
                            strainNames = new ArrayList<String>(phylo.getAllStrains());
                            Collections.sort(strainNames);
                            
                            outHeader = new String[strainNames.size() + 1];
                            for(int i = 0; i < strainNames.size(); i++)
                            {
                                outHeader[i] = strainNames.get(i);
                            }
                            outHeader[outHeader.length - 1] = "genomicIntervals";
                        }
                        
                        Set<BitSet> phyloSdps = phylo.sdps(strainNames, minMinorStrainCount);
                        for(BitSet sdp: phyloSdps)
                        {
                            List<GenomeInterval> sdpIntervals = sdpMap.get(sdp);
                            if(sdpIntervals == null)
                            {
                                sdpIntervals = new ArrayList<GenomeInterval>();
                                sdpMap.put(sdp, sdpIntervals);
                            }
                            sdpIntervals.add(genoInt);
                        }
                    }
                    
                    if(strainNames == null)
                    {
                        System.out.println(
                                "The input file " + inFileName +
                                " appears empty");
                        System.exit(-1);
                    }
                    else
                    {
                        // write the header row
                        FlatFileWriter ffw = new FlatFileWriter(
                                //new OutputStreamWriter(System.out),
                                new FileWriter(outFileName),
                                CommonFlatFileFormat.CSV_UNIX);
                        ffw.writeRow(outHeader);
                        int strainCount = strainNames.size();
                        String[] currOutRow = new String[strainNames.size() + 1];
                        for(Map.Entry<BitSet, List<GenomeInterval>> entry : sdpMap.entrySet())
                        {
                            BitSet currSDP = entry.getKey();
                            for(int i = 0; i < strainCount; i++)
                            {
                                currOutRow[i] = currSDP.get(i) ? "1" : "0";
                            }
                            
                            StringBuilder sb = new StringBuilder();
                            List<GenomeInterval> intervals = entry.getValue();
                            int intervalCount = intervals.size();
                            for(int i = 0; i < intervalCount; i++)
                            {
                                GenomeInterval interval = intervals.get(i);
                                sb.append(interval.getChrId());
                                sb.append(';');
                                sb.append(interval.getStartPositionBp());
                                sb.append(';');
                                sb.append(interval.getEndPositionBp());
                                if(i < intervalCount - 1)
                                {
                                    sb.append('|');
                                }
                            }
                            currOutRow[strainCount] = sb.toString();
                            ffw.writeRow(currOutRow);
                        }
                        
                        ffw.flush();
                    }
                }
            }
        }
        catch(ParseException ex)
        {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp("phylosdp", options);
            
            System.exit(-1);
        }
    }
}
