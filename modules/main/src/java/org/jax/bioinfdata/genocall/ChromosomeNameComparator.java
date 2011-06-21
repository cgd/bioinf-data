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

import java.io.Serializable;
import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A comparator for ordering chromosomes correctly according to their names.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ChromosomeNameComparator implements Comparator<String>, Serializable
{
    // after matching this pattern against a valid chromosome ID group 2
    // should contain the chromosome number (or in the case or X, Y and M
    // the chromosome letter)
    private static final Pattern CHR_PATTERN = Pattern.compile(
            "^(chromosome|chr)?\\s*(\\S+)$", Pattern.CASE_INSENSITIVE);
    
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -2436657796544414170L;

    /**
     * {@inheritDoc}
     */
    public int compare(String chrName1, String chrName2)
    {
        return chrToInt(chrName1) - chrToInt(chrName2);
    }
    
    private static int chrToInt(String chr)
    {
        Matcher matcher = CHR_PATTERN.matcher(chr);
        if(matcher.matches())
        {
            String matchedChr = matcher.group(2);
            // numbers come first followed by X, Y and M
            try
            {
                return Integer.parseInt(matchedChr);
            }
            catch(NumberFormatException ex)
            {
                matchedChr = matchedChr.toUpperCase();
                if(matchedChr.equals("X"))
                {
                    return Integer.MAX_VALUE - 2;
                }
                else if(matchedChr.equals("Y"))
                {
                    return Integer.MAX_VALUE - 1;
                }
                else if(matchedChr.equals("M"))
                {
                    return Integer.MAX_VALUE;
                }
                else
                {
                    throw new IllegalArgumentException(
                            "\"" + chr + "\" is not a valid chromosome name");
                }
            }
        }
        else
        {
            throw new IllegalArgumentException(
                    "\"" + chr + "\" is not a valid chromosome name");
        }
    }
}
