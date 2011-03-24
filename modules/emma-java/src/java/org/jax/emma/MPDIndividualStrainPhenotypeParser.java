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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jax.util.io.CharacterDelimitedParser;

/**
 * For parsing phenotype data in the 
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class MPDIndividualStrainPhenotypeParser
{
    /**
     * our comma separated file parser
     */
    private final CharacterDelimitedParser tabDelimitedParser =
        new CharacterDelimitedParser('\t');
    
    private static final String STRAIN_NAME_COLUMN_HEADER = "strain";
    
    private static final String PHENOTYPE_NAME_COLUMN_HEADER = "varname";
    
    private static final String SEX_COLUMN_HEADER = "sex";
    
    private static final String PHENOTYPE_VALUE_HEADER = "value";
    
    /**
     * Get available strain names from the stream.
     * @param is
     *          the stream to parse
     * @return
     *          the available phenotypes
     * @throws IOException
     *          if we catch an exception from the stream
     */
    public Set<String> parseAvailableStrainNames(InputStream is) throws IOException
    {
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        
        String[] headerArray = null;
        
        int sexColumnIndex = -1;
        int phenotypeNameIndex = -1;
        int phenotypeValueIndex = -1;
        int strainNameIndex = -1;
        
        do
        {
            String headerString = reader.readLine();
            
            if(headerString == null)
            {
                throw new IOException("failed to read header");
            }
            else if(!headerString.startsWith("#") && !headerString.startsWith("//"))
            {
                headerArray = this.tabDelimitedParser.parseCharacterDelimitedLine(
                        headerString);
                List<String> headerList = Arrays.asList(headerArray);
                
                sexColumnIndex = headerList.indexOf(SEX_COLUMN_HEADER);
                phenotypeNameIndex = headerList.indexOf(
                        PHENOTYPE_NAME_COLUMN_HEADER);
                phenotypeValueIndex = headerList.indexOf(
                        PHENOTYPE_VALUE_HEADER);
                strainNameIndex = headerList.indexOf(
                        STRAIN_NAME_COLUMN_HEADER);
            }
        } while(headerArray == null || sexColumnIndex == -1 ||
                phenotypeNameIndex == -1 || phenotypeValueIndex == -1 ||
                strainNameIndex == -1);
        
        Set<String> strainNames = new HashSet<String>();
        String[] currLine;
        while((currLine = this.tabDelimitedParser.parseCharacterDelimitedLine(reader)) != null)
        {
            strainNames.add(currLine[strainNameIndex]);
        }
        
        return strainNames;
    }
    
    /**
     * Get a set of all available phenotypes
     * @param is
     *          the input stream to parse
     * @return
     *          the available phenotypes
     * @throws IOException
     *          if we catch an exception
     */
    public Set<String> parseAvailablePhenotypes(InputStream is) throws IOException
    {
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        
        String[] headerArray = null;
        
        int sexColumnIndex = -1;
        int phenotypeNameIndex = -1;
        int phenotypeValueIndex = -1;
        int strainNameIndex = -1;
        
        do
        {
            String headerString = reader.readLine();
            
            if(headerString == null)
            {
                throw new IOException("failed to read header");
            }
            else if(!headerString.startsWith("#") && !headerString.startsWith("//"))
            {
                headerArray = this.tabDelimitedParser.parseCharacterDelimitedLine(
                        headerString);
                List<String> headerList = Arrays.asList(headerArray);
                
                sexColumnIndex = headerList.indexOf(SEX_COLUMN_HEADER);
                phenotypeNameIndex = headerList.indexOf(
                        PHENOTYPE_NAME_COLUMN_HEADER);
                phenotypeValueIndex = headerList.indexOf(
                        PHENOTYPE_VALUE_HEADER);
                strainNameIndex = headerList.indexOf(
                        STRAIN_NAME_COLUMN_HEADER);
            }
        } while(headerArray == null || sexColumnIndex == -1 ||
                phenotypeNameIndex == -1 || phenotypeValueIndex == -1 ||
                strainNameIndex == -1);
        
        Set<String> phenotypeNames = new HashSet<String>();
        String[] currLine;
        while((currLine = this.tabDelimitedParser.parseCharacterDelimitedLine(reader)) != null)
        {
            phenotypeNames.add(currLine[phenotypeNameIndex]);
        }
        
        return phenotypeNames;
    }
    
    /**
     * Parse phenotype data from the given stream.
     * @param phenotypeToParse
     *          the phenotype that we're reading in
     * @param inputStream
     *          the input stream
     * @param sexFilter
     *          the sex of the data we're parsing
     * @return
     *          the mapping of strains to phenotype values
     * @throws IOException
     *          if we catch an exception from the stream
     */
    public Map<String, List<Double>> parsePhenotypesFromStream(
            String phenotypeToParse,
            InputStream inputStream,
            SexFilter sexFilter)
    throws
            IOException
    {
        return this.parsePhenotypesFromStream(
                phenotypeToParse,
                inputStream,
                sexFilter,
                null);
    }
    
    /**
     * Parse phenotype data from the given stream
     * @param phenotypeToParse
     *          the phenotype to read
     * @param inputStream
     *          the input stream to read from
     * @param sexFilter
     *          the sex type that we're accepting
     * @param strainsToParse
     *          the strains that we should parse. to parse all strains use
     *          {@link #parsePhenotypesFromStream(String, InputStream, SexFilter)}
     *          instead
     * @return
     *          a mapping of strains to phenotype values
     * @throws IOException
     *          if we catch an exception from the stream
     */
    public Map<String, List<Double>> parsePhenotypesFromStream(
            String phenotypeToParse,
            InputStream inputStream,
            SexFilter sexFilter,
            Set<String> strainsToParse)
    throws
            IOException
    {
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        
        String[] headerArray = null;
        
        int sexColumnIndex = -1;
        int phenotypeNameIndex = -1;
        int phenotypeValueIndex = -1;
        int strainNameIndex = -1;
        
        do
        {
            String headerString = reader.readLine();
            
            if(headerString == null)
            {
                throw new IOException("failed to read header");
            }
            else if(!headerString.startsWith("#") && !headerString.startsWith("//"))
            {
                headerArray = this.tabDelimitedParser.parseCharacterDelimitedLine(
                        headerString);
                List<String> headerList = Arrays.asList(headerArray);
                
                sexColumnIndex = headerList.indexOf(SEX_COLUMN_HEADER);
                phenotypeNameIndex = headerList.indexOf(
                        PHENOTYPE_NAME_COLUMN_HEADER);
                phenotypeValueIndex = headerList.indexOf(
                        PHENOTYPE_VALUE_HEADER);
                strainNameIndex = headerList.indexOf(
                        STRAIN_NAME_COLUMN_HEADER);
            }
        } while(headerArray == null || sexColumnIndex == -1 ||
                phenotypeNameIndex == -1 || phenotypeValueIndex == -1 ||
                strainNameIndex == -1);
        
        Map<String, List<Double>> strainNameToPhenotypeValuesMap =
            new HashMap<String, List<Double>>();
        String[] currLine = new String[headerArray.length];
        while((currLine = this.tabDelimitedParser.parseCharacterDelimitedLine(reader)) != null)
        {
            // filter on values
            String sex = currLine[sexColumnIndex];
            if(sexFilter == SexFilter.ALLOW_FEMALE)
            {
                if(!sex.toLowerCase().startsWith("f"))
                {
                    continue;
                }
            }
            else if(sexFilter == SexFilter.ALLOW_MALE)
            {
                if(!sex.toLowerCase().startsWith("m"))
                {
                    continue;
                }
            }
            
            // get the phenotype name
            String phenotypeName = currLine[phenotypeNameIndex];
            if(!phenotypeName.equals(phenotypeToParse))
            {
                continue;
            }
            
            // get the strain name
            String currStrainName = currLine[strainNameIndex];
            if(strainsToParse != null && !strainsToParse.contains(currStrainName))
            {
                continue;
            }
            
            // add the phenotype value
            List<Double> phenotypesForCurrentStrain =
                strainNameToPhenotypeValuesMap.get(currStrainName);
            if(phenotypesForCurrentStrain == null)
            {
                phenotypesForCurrentStrain = new ArrayList<Double>();
                strainNameToPhenotypeValuesMap.put(
                        currStrainName,
                        phenotypesForCurrentStrain);
            }
            phenotypesForCurrentStrain.add(
                    Double.parseDouble(currLine[phenotypeValueIndex]));
        }
        
        return strainNameToPhenotypeValuesMap;
    }
}
