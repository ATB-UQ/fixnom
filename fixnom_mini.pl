#!/usr/bin/perl
#
#Program to rename MOL2 file according to IUPAC
#
use Getopt::Long;

#READ OPTIONS ****************************
GetOptions('ca=s' => \$ca,
           'n=s' => \$n,
           'c=s' => \$c,
           'i=s' => \$inputfile,
           'o=s' => \$outputfile
			);
# -c C6 -x1 C7 -x2 N2
#READ LIB FILE **************************
open(LIB, $inputfile) || die "Error: Could not open input library file.\nUsage\n ./fixnom.pl -i input.lib -o output.lib -ca C3 -c C4 -n N1\n ";		# Open the input file
@rlib = <LIB>;										# Read it into an array
close(LIB);
# ***********************
# READ GREEK NAMES
my @gnam  =('A', 'B', 'G', 'D', 'E', 'Z', 'H', 'T', 'I', 'K', 'L', 'M', 'N', 'X', 'O', 'P', 'R', 'S', 'J', 'U', 'F', 'C', 'Y', 'W');
$gind=0;    # index to keep track of names
#********


#********

$foundres=$found_n=$found_c=$found_o=0;

	$j=$l=0;
    $i=1;
    $x=0;
$atstart = $hlin = $l = $k = 0;
 $nterm=$cterm=0;
$at_n=$at_c=$at_o=-1;
foreach $line (@rlib)
	{
		@fields = split (' ', $line);
        chomp($line);
        $temp_crn[$i] = 0;
        
# Flagging we are in a residue
       # ATOMTYPES     1      0.10
       # 1 CHLOR      1.00    1   17
        if ($fields[0] =~ '@<TRIPOS>MOLECULE')
        {
            $inflin=$i+2;
        }
        if ($fields[0] =~ '@<TRIPOS>ATOM')
        {
            $inat=$i+1;
            $endat=$i+$numat;
        }
        if ($fields[0] =~ '@<TRIPOS>BOND')
        {
            $inbnd=$i+1;
            $endbnd=$i+$numbnd;
        }
        if ($i == $inflin)
        {
            $numat = $fields[0];
            $numbnd = $fields[1];
            #printf "%i %i\n", $numat, $numbnd;
        }
#        unless ($fields[5] =~ 'H' && length($fields[5]=1))
#        {
		 if ($i >= $inat && $i <= $endat )
		 {
            #printf "%s\n", $line;
            $l++;
            $temp_c1[$l] = $fields[0];
            $temp_c2[$l] = $fields[1];
            @fields2 = split (/\./, $fields[5]);
            $temp_c3[$l] = $fields2[0];
            $temp_c4[$l] = 0;
            $temp_c6[$l] = 0;
            $temp_c5[$l]= $l; # pointer list
            $fil3_c[$l]=$line;
            $fil2_c[$i]=3;

            #printf "%s %s %s %s\n", $fields[0], $fields[1], $fields[5], $fields2[0];

            if ($ca =~ $temp_c2[$l]  &&  length($ca) == length($temp_c2[$l]))
            {
                $at_ca=$temp_c1[$l];
                #printf "%s %s %s %s\n", $fields[0], $fields[1], $fields[5], $ca;
            }
            if ($c =~ $temp_c2[$l]  &&  length($c) == length($temp_c2[$l]))
            {
                $at_c=$temp_c1[$l];
                #printf "%s %s %s %s\n", $fields[0], $fields[1], $fields[5], $ca;
            }
            if ($n =~ $temp_c2[$l]  &&  length($n) == length($temp_c2[$l]))
            {
                $at_n=$temp_c1[$l];
                #printf "%s %s %s %s\n", $fields[0], $fields[1], $fields[5], $ca;
            }
         }elsif ($i >= $inbnd && $i <= $endbnd)
         {
            $k++;
            $temp_b1[$k] = $fields[1];
            $temp_b2[$k] = $fields[2];
            $temp_b3[$k] = $fields[3];
            $fil2_c[$i]=2;

         }else{
        #printf "%s\n", $line;
             $fil2_c[$i]=1;
        }
        $fil_c[$i]=$line;           #stores each line in file to array $fil_c
        $i++;                       #counts rows - after row $fields[2] the atoms starts
	}   # end of reading file

$numlin		= $i-1;
$numhat     = $l;

#################################get NH AND CO ########################
for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c1[$i] == $at_c)
    {
        for ($j = 1 ; $j <= $numbnd ; $j++)
        {
            unless ($temp_b1[$j] == $at_ca || $temp_b2[$j] == $at_ca)
            {
                if($temp_b1[$j] == $temp_c1[$i] && $temp_c3[$temp_b2[$j]] =~ 'O')
                {
                    $at_co = $temp_b2[$j] ;
                }
                if($temp_b2[$j] == $temp_c1[$i] && $temp_c3[$temp_b1[$j]] =~ 'O')
                {
                    $at_co = $temp_b1[$j] ;
                }
            }
        }
    }
    if($temp_c1[$i] == $at_n)
    {
        for ($j = 1 ; $j <= $numbnd ; $j++)
        {
            unless ($temp_b1[$j] == $at_ca || $temp_b2[$j] == $at_ca)
            {
                if($temp_b1[$j] == $temp_c1[$i] && $temp_c3[$temp_b2[$j]] =~ 'H')
                {
                    $at_hn = $temp_b2[$j] ;
                }
                if($temp_b2[$j] == $temp_c1[$i] && $temp_c3[$temp_b1[$j]] =~ 'H')
                {
                    $at_hn = $temp_b1[$j] ;
                }
            }
        }
    }
}
################################# CUT OUT CO SIDE OF THINGS ######################
# NEED TO reduce matrix further to cut out N/C connections.
# assumes index goes form 1 to end.
# Sort so that CO is on top
for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c1[$i] == $at_c && $temp_c1[1] != $temp_c1[$i])
    {
        ($temp_c1[1], $temp_c1[$i]) = ($temp_c1[$i], $temp_c1[1]);
        ($temp_c2[1], $temp_c2[$i]) = ($temp_c2[$i], $temp_c2[1]);
        ($temp_c3[1], $temp_c3[$i]) = ($temp_c3[$i], $temp_c3[1]);
        ($temp_c4[1], $temp_c4[$i]) = ($temp_c4[$i], $temp_c4[1]);

    }
}
#printf "oxygen %i\n", $at_co;
$temp_c4[1]=1;
$t1=1;
for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c4[$i] == 1)
    {
        for ($j = 1 ; $j <= $numbnd ; $j++)
        {
#            printf "%i %i %i\n", $temp_c1[$i], $temp_b1[$j], $at_ca;

            unless ($temp_b1[$j] == $at_ca || $temp_b2[$j] == $at_ca)
            {
                if($temp_b1[$j] == $temp_c1[$i] )
                {
                    #printf "%i %s %i\n", $temp_b2[$j], $temp_c2[$temp_b2[$j]], $i;
                    #need to match b2 to atom list
                    for ($i2=$i+1 ; $i2 <= $numhat ; $i2++)      ### just for printing
                    {
                        if ($temp_c4[$i2] == 0 && $temp_c1[$i2] == $temp_b2[$j])
                        {
                            $temp_c4[$i2]=1;
                        }
                    }
                        
                    #$temp_c4[$temp_b2[$j]]=$t1;
                    #$t1=1;
                    #$a1p=$i;
                }
                if($temp_b2[$j] == $temp_c1[$i] )
                {
                    #$a2p=$i;
                    #printf "1 %i %s %i\n", $temp_b1[$j], $temp_c2[$temp_b1[$j]], $i;
                    #$temp_c4[$temp_b1[$j]]=$t1;
                    for ($i2=$i+1 ; $i2 <= $numhat ; $i2++)      ### just for printing
                    {
                        if ($temp_c4[$i2] == 0 && $temp_c1[$i2] == $temp_b1[$j])
                        {
                            $temp_c4[$i2]=1;
                        }
                    }
                }
                
            }
        }
    }
    ## sort the array but record the ordered indexes.
    for ($j=$i+1 ; $j <= $numhat ; $j++)      ### just for printing
    {
        for ($k=$j+1 ; $k < $numhat ; $k++)      ### just for printing
        {
            if($temp_c4[$k] > $temp_c4[$j])
            {
                ($temp_c1[$j], $temp_c1[$k]) = ($temp_c1[$k], $temp_c1[$j]);
                ($temp_c2[$j], $temp_c2[$k]) = ($temp_c2[$k], $temp_c2[$j]);
                ($temp_c3[$j], $temp_c3[$k]) = ($temp_c3[$k], $temp_c3[$j]);
                ($temp_c4[$j], $temp_c4[$k]) = ($temp_c4[$k], $temp_c4[$j]);
            }
    
        }
    }
#     print "@$_ \n" for \( @temp_c1, @temp_c2, @temp_c4 );
}

################################# CUT OUT CO SIDE OF THINGS ######################

################################# CUT OUT N SIDE OF THINGS ######################
# NEED TO reduce matrix further to cut out N/C connections.
# assumes index goes form 1 to end.
# Sort so that CO is on top
for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c1[$i] == $at_n && $temp_c1[1] != $temp_c1[$i])
    {
        ($temp_c1[1], $temp_c1[$i]) = ($temp_c1[$i], $temp_c1[1]);
        ($temp_c2[1], $temp_c2[$i]) = ($temp_c2[$i], $temp_c2[1]);
        ($temp_c3[1], $temp_c3[$i]) = ($temp_c3[$i], $temp_c3[1]);
        ($temp_c4[1], $temp_c4[$i]) = ($temp_c4[$i], $temp_c4[1]);

    }
}
#printf "oxygen %i\n", $at_co;
$temp_c4[1]=-1;
for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c4[$i] == -1)
    {
        for ($j = 1 ; $j <= $numbnd ; $j++)
        {
#            printf "%i %i %i\n", $temp_c1[$i], $temp_b1[$j], $at_ca;

            unless ($temp_b1[$j] == $at_ca || $temp_b2[$j] == $at_ca)
            {
                if($temp_b1[$j] == $temp_c1[$i] )
                {
                    #printf "%i %s %i\n", $temp_b2[$j], $temp_c2[$temp_b2[$j]], $i;
                    #need to match b2 to atom list
                    for ($i2=$i+1 ; $i2 <= $numhat ; $i2++)      ### just for printing
                    {
                        if ($temp_c4[$i2] == 0 && $temp_c1[$i2] == $temp_b2[$j])
                        {
                            $temp_c4[$i2]=-1;
                        }
                    }
                        
                }
                if($temp_b2[$j] == $temp_c1[$i] )
                {

                    for ($i2=$i+1 ; $i2 <= $numhat ; $i2++)      ### just for printing
                    {
                        if ($temp_c4[$i2] == 0 && $temp_c1[$i2] == $temp_b1[$j])
                        {
                            $temp_c4[$i2]=-1;
                        }
                    }
                }
                
            }
        }
    }
    ## sort the array but record the ordered indexes.
    for ($j=$i+1 ; $j <= $numhat ; $j++)      ### just for printing
    {
        for ($k=$j+1 ; $k < $numhat ; $k++)      ### just for printing
        {
            if($temp_c4[$k] < $temp_c4[$j])
            {
                ($temp_c1[$j], $temp_c1[$k]) = ($temp_c1[$k], $temp_c1[$j]);
                ($temp_c2[$j], $temp_c2[$k]) = ($temp_c2[$k], $temp_c2[$j]);
                ($temp_c3[$j], $temp_c3[$k]) = ($temp_c3[$k], $temp_c3[$j]);
                ($temp_c4[$j], $temp_c4[$k]) = ($temp_c4[$k], $temp_c4[$j]);
            }
    
        }
    }
#     print "@$_ \n" for \( @temp_c1, @temp_c2, @temp_c4 );
}


$temp_c1[0]=-9999;

my @order = sort{ $temp_c1[ $a ] <=> $temp_c1[ $b ] } @0 .. $#temp_c1;

@temp_c1 = @temp_c1[ @order ];
@temp_c2 = @temp_c2[ @order ];
@temp_c3 = @temp_c3[ @order ];
@temp_c4 = @temp_c4[ @order ];



#print "@$_ \n" for \( @temp_c1, @temp_c2, @temp_c4 );

#$fil3_c[$l]

printf "N-term block.\n";
for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c4[$i] == -1 && $temp_c1[$i] != $at_n && $temp_c1[$i] != $at_hn)
    {
        printf "%i %s\n", $temp_c1[$i], $temp_c2[$i];
        $temp_c4[$i] =-1;
    }
}
printf "Central block.\n";
$l=0;
for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c4[$i] == 0 || $temp_c1[$i] == $at_c || $temp_c1[$i] == $at_co || $temp_c1[$i] == $at_hn || $temp_c1[$i] == $at_n)
    {
        printf "%i %s\n", $temp_c1[$i], $temp_c2[$i];
        unless ($temp_c3[$i] =~ 'H')
        {        if($temp_c4[$i] == 0 )
            {
                $l++;
                $lmat[1][$l]=$temp_c1[$i];
                $lmat[12][$l]=$temp_c2[$i];
                $lmat[2][$l]=$temp_c3[$i];
            }
        }
        $temp_c4[$i]=0;     # to bring in the NH and CO into central block
    }

}
printf "C-term block.\n";
for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    
    if($temp_c4[$i] == 1 && $temp_c1[$i] != $at_c && $temp_c1[$i] != $at_co)
    {
        printf "%i %s\n", $temp_c1[$i], $temp_c2[$i];
        $temp_c4[$i] =1;
    }
}
$onumhat=$numhat;
$numhat=$l;
################################# CUT OUT N SIDE OF THINGS ######################

############### PRINT #############


$temp_c4[0]=-9999;

my @order = sort{ $temp_c4[ $a ] <=> $temp_c4[ $b ] } @0 .. $#temp_c4;

@temp_c1 = @temp_c1[ @order ];
@temp_c2 = @temp_c2[ @order ];
@temp_c3 = @temp_c3[ @order ];
@temp_c4 = @temp_c4[ @order ];
@fil3_c  = @fil3_c[ @order ];
$mol2 = sprintf("block_%s", $inputfile);

open(OUT2, '>', $mol2) || die "Can't open output file";
#        printf OUT "%s",$head[$i];

$endb=0;
for ($j = 1 ; $j < $numlin ;  $j++)
{
     if ($fil2_c[$j] == 3)
     {
        if($endb == 0)
        {
            for ($i=1 ; $i <= $onumhat ; $i++)      ### just for printing
            {
                printf OUT2 "%s\n", $fil3_c[$i];
            }
        }
     $endb=1;
     }else{
     printf OUT2 "%s\n", $fil_c[$j];
     }
}
close OUT2;
     #$fil3_c[$l]

$temp_c1[0]=-9999;

my @order = sort{ $temp_c1[ $a ] <=> $temp_c1[ $b ] } @0 .. $#temp_c1;

@temp_c1 = @temp_c1[ @order ];
@temp_c2 = @temp_c2[ @order ];
@temp_c3 = @temp_c3[ @order ];
@temp_c4 = @temp_c4[ @order ];




############### PRINT #############



### REDUCED MATRIX FROM mol2 FILE CREATED ########
#exit;



for ($i = 1 ; $i <= $numbnd ; $i++)
{
    $a1 = $temp_b1[$i];
    $a2 = $temp_b2[$i];
    $t1=$t2=0;
    for ($j=1 ; $j <= $numhat ; $j++)      ### just for printing
    {
        if($a1 == $lmat[1][$j])
        {
            $t1=1;
            $a1p=$j;
        }
        if($a2 == $lmat[1][$j])
        {
            $t2=1;
            $a2p=$j;
        }
    }
    if ($t1 == $t2 && $t1 == 1)
    {
      #  printf "bond between two hetero %i %i %i %i\n", $a1, $a2, $a1p, $a2p;
        $bond_lst_l[1][$a1p]++;
        $bond_lst_l[1][$a2p]++;
        $l1=$bond_lst_l[1][$a1p];
        $l2=$bond_lst_l[1][$a2p];
        $master_c[$a1p]=$l1;
        $master_c[$a2p]=$l2;

        $bond_list[1][$a1p][$l1] = $a2;
        $bond_list[1][$a2p][$l2] = $a1;
        $master_l[$a1p][$l1] = $a2;   # first entry is itself
        $master_l[$a2p][$l2] = $a1;   # first entry is itself

        $bond_list_ptr[1][$a1p][$l1]=$a2p;
        $bond_list_ptr[1][$a2p][$l2]=$a1p;
        
    }
}

for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{

    for ($x = 1; $x < 2 ; $x++)     # 1 - bond matrix check
    {
        for ($j=1 ; $j <= $bond_lst_l[$x][$i] ; $j++)       # how many entries in x-bond list - go through each
        {
            #printf "%s id = %i, i = %i, j= %i, x = %i, list:  %i\n", $lmat[12][$i], $lmat[1][$i], $i, $j, $x, $bond_list[$x][$i][$j];
        }
    }
#    printf "%i %i %s\n", $lmat[1][$i], $bond_lst_l[1][$i], $lmat[2][$i];
#    $fpr_i[$i]=0;
}
#exit;


if ($at_ca < 1 || $at_n < 1 )
{
    
   printf "Error – backbone C alpha and/or N not found. \n";
    exit;

}
$nhat=$numhat;
#    printf "mid %i %i %i %i\n", $at_ca, $at_c, $nhat,  $lmat[1][10];


# condensed array created
$nhat = $nhat+1; # cut out carbonyl carbon start at 1 to cut out backbone N.
# if we know what at_ca is we can skip to here.
#printf "ats %i %i %i %i\n",  $nhat, $numhat, $at_ca, $lmat[1][1];
for ($i=1 ; $i < $nhat ; $i++)
{
    if ($lmat[1][$i] == $at_ca)
    {
        for ($j=1; $j <= 12 ; $j++)      # c-alpha isn't the first atom reorder to make this the case
        {
            $tvec[$j] = $lmat[$j][1];
            $lmat[$j][1] = $lmat[$j][$i];
            $lmat[$j][$i] = $tvec[$j];
        }
        
        for ($j=1 ; $j <= $bond_lst_l[1][1] ; $j++)       # how many entries in x-bond list - go through each
        {
                 $tvec[$j]=    $bond_list[1][1][$j];
                $tvecp[$j]=$bond_list_ptr[1][1][$j];
        }
        $temp_l=$bond_lst_l[1][1];

        $bond_lst_l[1][1]=$bond_lst_l[1][$i];
        for ($j=1 ; $j <= $bond_lst_l[1][1] ; $j++)       # how many entries in x-bond list - go through each
        {
                $bond_list[1][1][$j]=    $bond_list[1][$i][$j];
            $bond_list_ptr[1][1][$j]=$bond_list_ptr[1][$i][$j];
        }

        $bond_lst_l[1][$i]=$temp_l;
        for ($j=1 ; $j <= $bond_lst_l[1][$i] ; $j++)       # how many entries in x-bond list - go through each
        {
                $bond_list[1][$i][$j]= $tvec[$j];
            $bond_list_ptr[1][$i][$j]=$tvecp[$j];
        }
    }
    
}


#exit;
# initialise matrix - matrix goes from N amide residue onwards until carbonyl carbon and then stops - carbony carbon must be last residue - can be cut here.

for ($i=1 ; $i < $nhat ; $i++)
{
    $attyp = $lmat[2][$i];          # get atom number  - used for priority
    &getatno($attyp, $atno);
    $lmat[10][$i]=$atno;
}
#   FORM CHAINS
#   CREATE PAIRWISE CONNECTIONS GOING to CA.
$nbonds=0;
$b=1;
$bc=0;
$fs=0;
for ($k=1 ; $k < $nhat-1 ; $k++)
{
    #$nbonds=$nbonds+1;
    $ptr=$lmat[1][$k];  #   [sets pointer to c-alpha position at first round]
    $lmat[7][$k]=1;     #   Atom is being considered/has been considered
    $b=0;
    for ($i=$k+1 ; $i < $nhat ; $i++)
        {

            if ($lmat[7][$i] == 1) {next} ;     # if atom has been seen previously skip to next element in for loop

            for ($j=1 ; $j <= $bond_lst_l[1][$i] ; $j++)       # how many entries in x-bond list - go through each
            {
                if($bond_list[1][$i][$j] == $ptr)       # something is connected to this point along the chain
                {
                    $lmat[9][$i]=$ptr;          #atom number of hit
                    $lmat[7][$i]=1;
                    last;                      # removed last so that priority (atom number can be carried forward)
                }
            }
        }
########  START SORTING       need to sort based on $lmat[7] from position $k +1        re-order based on distance to CA - i.e. seen
    for ($si=$k+1 ; $si < $nhat-1 ; $si++)
        {
            for ($sj = $si+1 ; $sj < $nhat ; $sj++)
            {
                if ($lmat[7][$si] < $lmat[7][$sj])
                {
                    for ($l=1; $l <= 12 ; $l++)
                    {
                                $tvec[$l]           = $lmat[$l][$si];
                                $lmat[$l][$si]      = $lmat[$l][$sj];
                                $lmat[$l][$sj]       = $tvec[$l];
                    }
                    
                    for ($j=1 ; $j <= $bond_lst_l[1][$si] ; $j++)       # how many entries in x-bond list - go through each
                    {
                         $tvec[$j] =    $bond_list[1][$si][$j];
                    }
                    $temp_l=$bond_lst_l[1][$si];

                    $bond_lst_l[1][$si]=$bond_lst_l[1][$sj];
                    for ($j=1 ; $j <= $bond_lst_l[1][$si] ; $j++)       # how many entries in x-bond list - go through each
                    {
                            $bond_list[1][$si][$j]=    $bond_list[1][$sj][$j];
                    }

                    $bond_lst_l[1][$sj]=$temp_l;
                    for ($j=1 ; $j <= $bond_lst_l[1][$sj] ; $j++)       # how many entries in x-bond list - go through each
                    {
                            $bond_list[1][$sj][$j]= $tvec[$j];
                    }
                }
                
            }
        }
#########END SORTING ROUTINE
}


for ($i=1 ; $i < $nhat ; $i++)      ### UPDATE PTR LIST
{

    for ($x = 1; $x < $nhat ; $x++)
    {
        for ($j=1 ; $j <= $bond_lst_l[$x][$i] ; $j++)       # how many entries in x-bond list - go through each
        {
            for ($i2=1 ; $i2 < $nhat ; $i2++)      ### just for printing
            {
                if ($bond_list[$x][$i][$j] == $lmat[1][$i2])
                {
                    $bond_list_ptr[$x][$i][$j] = $i2;
                }
            }
            #printf "NEWxx id = %i, i = %i, j= %i, x = %i, list:  %i ptr %i\n", $lmat[1][$i], $i, $j, $x, $bond_list[$x][$i][$j], $bond_list_ptr[$x][$i][$j];
        }
    }
#    printf "%i %i %s\n", $lmat[1][$i], $bond_lst_l[1][$i], $lmat[2][$i];
#    $fpr_i[$i]=0;
}#

#exit;

################## FOR CIP rule need to add a bunch of carbons now. These will need to be checked but not priorities or chained themselves.
### THIS IS NOT WORKING AND Therefore $dum = 0.
$dum=0;
#
#
$test=1;
if ($test == 1)
{
    
    for ($i = 1 ; $i <= $numbnd ; $i++)
    {
        if ($temp_b3[$i] =~ 'ar')
        {
            #printf "%i %i %s %i %i\n", $temp_b1[$i], $temp_b2[$i], $temp_b3[$i], $temp_c6[$temp_b1[$i]], $temp_c6[$temp_b2[$i]] ;

            $temp_b3[$i]=2;
            $temp_c6[$temp_b1[$i]]++;
            $temp_c6[$temp_b2[$i]]++;       # contineu from here
            

        }
        if ($temp_b3[$i] > 1)
        {
            #printf "%i %i\n", $temp_b1[$i], $temp_b3[$i];
            # find the two atoms and add in the dummy for each
                #printf "I'm here %i with %i %i\n", $lmat[1][$j], $temp_b2[$i], $temp_b1[$i]  ;

                    #printf "B1 I'm here %i with %i\n", $lmat[1][$j], $temp_b2[$i] ;
                    $a1p=$temp_b1[$i];
                    $a2p=$temp_b2[$i];

           # printf "ATOM %i IS A %s CONNECTED TO ATOM %i IS A %s WITH A %i BONDS \n", $temp_c1[$a1p], $temp_c3[$a1p], $temp_c1[$a2p], $temp_c3[$a2p], $temp_b3[$i] ;
           # printf "ATOM %i IS A %s CONNECTED TO ATOM %i IS A %s WITH A %i BONDS \n", $temp_c1[$a2p], $temp_c3[$a2p], $temp_c1[$a1p], $temp_c3[$a1p], $temp_b3[$i] ;

            $t1=$t2=0;
            for ($i2 = 1 ; $i2 < $nhat ; $i2++)
            {
                if($temp_c1[$a1p]==$lmat[1][$i2])
                {
                    $id1 = $i2;      #
                    $t1=1;
                }
                if($temp_c1[$a2p]==$lmat[1][$i2])
                {
                    $id2 = $i2;      #
                    $t2=1;
                }
            }
 
            if ($t1 == $t2 && $t1 == 1)
            {
             if ($temp_c6[$temp_b1[$i]] <= 1)
             {
                 for ($l = 1 ; $l < $temp_b3[$i] ; $l++)
                 {

                        $lmat[1][$nhat+$dum]=1000+$dum;         #atom ID
                        $lmat[2][$nhat+$dum]=$lmat[2][$id2];     #atom TYPE
                        $lmat[3][$nhat+$dum]=$lmat[1][$id1];      #atom CONNECTION
                        $lmat[4][$nhat+$dum]=0;                 #atom C2
                        $lmat[5][$nhat+$dum]=0;                 #atom C3
                        $lmat[6][$nhat+$dum]=0;                 #atom C4
                        $lmat[7][$nhat+$dum]=$lmat[7][$id1];      #atom CHAIN ID
                        $lmat[8][$nhat+$dum]=0;                 #atom C4
                        $lmat[9][$nhat+$dum]=$lmat[1][$id1];      #PAIRWISE CONNECTION (BACK TO CA)
                        $lmat[10][$nhat+$dum]=$lmat[10][$id2];    #atom NUMBER - FOR PRIORITY
                     #printf "dummy ATOM %i IS A %s CONNECTED TO ATOM %i WITH THE %i BOND. ch id %i p %i\n", $nhat+$dum, $lmat[2][$nhat+$dum], $lmat[9][$nhat+$dum], $l, $lmat[3][$nhat+$dum], $lmat[10][$nhat+$dum];

                     $dum++;
                 }
             }
             if ($temp_c6[$temp_b2[$i]] <= 1)
             {
                 for ($l = 1 ; $l < $temp_b3[$i] ; $l++)
                 {

                           $lmat[1][$nhat+$dum]=1000+$dum;         #atom ID
                           $lmat[2][$nhat+$dum]=$lmat[2][$id1];     #atom TYPE
                           $lmat[3][$nhat+$dum]=$lmat[1][$id2];      #atom CONNECTION
                           $lmat[4][$nhat+$dum]=0;                 #atom C2
                           $lmat[5][$nhat+$dum]=0;                 #atom C3
                           $lmat[6][$nhat+$dum]=0;                 #atom C4
                           $lmat[7][$nhat+$dum]=$lmat[7][$id2];      #atom CHAIN ID
                           $lmat[8][$nhat+$dum]=0;                 #atom C4
                           $lmat[9][$nhat+$dum]=$lmat[1][$id2];      #PAIRWISE CONNECTION (BACK TO CA)
                           $lmat[10][$nhat+$dum]=$lmat[10][$id1];    #atom NUMBER - FOR PRIORITY
                    #printf "dummy2 ATOM %i IS A %s CONNECTED TO ATOM %i WITH THE %i BOND. ch id %i p %i\n", $nhat+$dum, $lmat[2][$nhat+$dum], $lmat[9][$nhat+$dum], $l, $lmat[3][$nhat+$dum], $lmat[10][$nhat+$dum];

                    $dum++;
                    }
                
                }
            }
            
        }

    }
    $nhat=$nhat+$dum;

}
    
#    exit;
 

$test=0;
$test=$nhat;

#####################################$$$$$$$$$$$$$$$$$$$$$$$$$$
#################  all 1 bond 2 bond 3 bond etc. connectivities for all atoms. - to be used later for assigning priority
for ($i=1 ; $i < $nhat ; $i++)
{
#    printf "%i atoms 1 bond away from %s: \n", $lmat[1][$i], $lmat[12][$i];
    $l=0;
    $master_l[$i][$l] = $lmat[1][$i];   # first entry is itself

        for ($j=1 ; $j <= $bond_lst_l[1][$i] ; $j++)       # how many entries in x-bond list - go through each
        {
#            printf "%i ", $bond_list[1][$i][$j];
            $master_l[$i][$j] = $bond_list[1][$i][$j];   # first entry is itself
        }
    
#    printf "\n";
    $master_c[$i]=$bond_lst_l[1][$i];
}
####   USING ABOVE 1 BOND ARRAY - CREATE GENERAL >2 ATOMS AWAY ARRAYS THAT ARE X ATOMS AWAY FROM EACH ENTRY
for ($i=1 ; $i < $nhat ; $i++)
{
    $l=$master_c[$i];
    for ($x = 1; $x < $nhat ; $x++)
    {
        $bond_lst_l[$x+1][$i] =0;
        for ($j=1 ; $j <= $bond_lst_l[$x][$i] ; $j++)       # how many entries in x-bond list - go through each
        {
            $ptr=$bond_list_ptr[$x][$i][$j];                # pointer to current point along list
            #printf "i = %i, j= %i, x = %i, list:  %i, ll: %i, ML: %i\n", $i, $j, $x, $bond_list[$x][$i][$j], $bond_lst_l[1][$ptr], $l;
#            $bond_lst_l[$x+1][$i] = $bond_lst_l[$x][$ptr]+$bond_lst_l[$x+1][$i] ;
                for ($j2=1 ; $j2 <= $bond_lst_l[1][$ptr] ; $j2++)       # how many entries in x-bond list - go through each
                {
                    $skip = 0;
                    for ($l2 = 0 ; $l2 <= $l ; $l2++)
                    {
                        if ($bond_list[1][$ptr][$j2] == $master_l[$i][$l2] || $bond_list[1][$ptr][$j2] == $lmat[1][$i])
                        {
                            $skip = 1;
                            last;
                        }
                        #printf "l2 = %i, l= %i, master = %i, skip: %i, j2: %i, ptr: %i\n", $l2, $l, $master_l[$i][$l2], $skip, $j2, $bond_list[1][$ptr][$j2];
                    }
                    ### loop through bondlist x and see if current entry is present
                    if ($skip == 0)
                    {
                        $l++;
                        $bond_lst_l[$x+1][$i]++;
                        $j3=$bond_lst_l[$x+1][$i];
                             $bond_list[$x+1][$i][$j3]=    $bond_list[1][$ptr][$j2];
                         $bond_list_ptr[$x+1][$i][$j3]=$bond_list_ptr[1][$ptr][$j2];
                        $master_l[$i][$l]=$bond_list[1][$ptr][$j2];       # contains all past entries
                    }
                }

        }
        if ($bond_lst_l[$x+1][$i] == 0)
        {
            last;
        }
    }
    $master_c[$i]=$l;
}
for ($i=1 ; $i < $nhat ; $i++)      ### just for printing
{
        $bond_list[0][$i][1]=$lmat[1][$i];
    $bond_list_ptr[0][$i][1]=$i;
       $bond_lst_l[0][$i]=1;
    for ($x = 1; $x < $nhat ; $x++)
    {
        for ($j=1 ; $j <= $bond_lst_l[$x][$i] ; $j++)       # how many entries in x-bond list - go through each
        {
#            printf "err id = %i, i = %i, j= %i, x = %i, list:  %i, ptr %i\n", $lmat[1][$i], $i, $j, $x, $bond_list[$x][$i][$j], $bond_list_ptr[$x][$i][$j];
        }
    }
    $fpr_i[$i]=0;
}
#exit;
####################### finish building tree! ^=##########
#####################################$$$$$$$$$$$$$$$$$$$$$$$$$$

##################### ASSIGN CHAIN ID
$chains=0;
$t_chains=1;
$c_chains=0; #current and total chains
$cxl[$chains]=0;
$lmat[8][1]=1;  # calpha is chain 1
$lmat[11][1]=0;
#       $lmat[8] will be used here to index chains
#$lmat[10][5]=8;
#printf "Modified: %i\n", $lmat[1][5];
#$lmat[10][10]=8;
#printf "DUM %i\n", $dum;
$nhat=$nhat-$dum;
for ($i=1 ; $i < $nhat-1 ; $i++)
{
#    $attyp = $lmat[2][$i];
#    &getatno($attyp, $atno);
#    $lmat[10][$i]=$atno;
    $t_chains=$t_chains+$c_chains;
    $c_chains=0;
    $found=0;
    $addchain=0;
    $fpr_i[$i]=0;

############        ASSIGNING PRIORITY  ############
                if ($bond_lst_l[1][$i]>2)               # if more than 2 attached heteroatmos, then this is a branch
                {
                    
                    for ($j2=1 ; $j2 <= $bond_lst_l[1][$i] ; $j2++)         # check how many branching atoms > $i - > if less than $i then not branch
                    {
                        #printf "PTR %i\n", $lmat[1][$bond_list_ptr[1][$i][$j2]];

                        if ($bond_list_ptr[1][$i][$j2] > $i && $lmat[1][$bond_list_ptr[1][$i][$j2]] < 1000)                # since things are sorted above if the id is less than current ID then it can be ignored
                        #if ($bond_list_ptr[1][$i][$j2] > $i )                # since things are sorted above if the id is less than current ID then it can be ignored
                        {
                            $addchain++;                                    # if more than 1 you really do have a new branch - addchain tells you how many chains needed
                            $fpr[$addchain]=$bond_list_ptr[1][$i][$j2];     # ptr of first atom in new chain
                            $cord[$addchain]=$addchain;                     # chain value/number
                                                                            # add chain tells you how many chains that need to be prioritised.
                        }
                    }
                }
#       i have list of atoms that need to be prioritised
        if($addchain > 1)
        {
################ Check priority of hits and sort based on hit #####################
###########      Check priority of atom itself.                 ######## needs to be done separately as it is not part of the bond_list.
            
            
            
            
            
####### NEED TO CREATE A MATRIX THAT LOOKS LIKE THIS
#           C1  C2  C3  ...C_ADDCHAIN
#   SHELL 0
#   SHELL 1
#   SHELL 2
#   .
#   .
#   ALL ENTRIES ARE 0
#   VALUES IN MATRIX ARE ATOM NUMBERS
#   for each of the atoms at position SHELL 0 write out the priorities.
#   for each of these atoms
            for ($p = 1 ; $p <= $addchain ; $p++)
            {
                $x=0;
                $i2 = $fpr[$p];
                $pr_mat[$p][$x] =  $lmat[10][$i2];
                for ($x = 1; $x < $nhat ; $x++)
                {
                    $pr_mat[$p][$x] =  0;
                    if ($bond_lst_l[$x][$i2] > 0)
                    {

                        $i3             = $bond_list_ptr[$x][$i2][1];
                        $maxp           = $lmat[10][$i3];                              # the first entry
                        $pr_mat[$p][$x] = $lmat[10][$i3];
                       # printf "Atom i: %i Bonds away: %i Priority: %i\n", $lmat[1][$fpr[$p]], $x, $lmat[10][$i3], ;

                        for ($j=1 ; $j <= $bond_lst_l[$x][$i2] ; $j++)       # how many entries in x-bond list - go through each
                        {
                            $i3=$bond_list_ptr[$x][$i2][$j];
                            if ($lmat[10][$i3] > $maxp )
                            {
                                $pr_mat[$p][$x] = $lmat[10][$i3];
                                $maxp           = $lmat[10][$i3];
                            }
                        }
                    }
                }
                $fpr_i[$i]=0;
            }
            for ($p = 1 ; $p <= $addchain-1 ; $p++)
            {
                for ($p2 = 1+$p ; $p2 <= $addchain ; $p2++)
                {
                    for ($x = 0; $x < $nhat ; $x++)
                    {
                        #printf "Atom i: %i Bonds away: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$p]], $x, $pr_mat[$p][$x], $p;
                        #printf "CF A i: %i Bonds away: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$p2]], $x, $pr_mat[$p2][$x], $p2;
                        #if ($pr_mat[$p][$x] == 0) {printf "should exit here saying that assignment is arbitrary\n"};
####################
#                        THIS IS WHERE THE GEOMETRIC CIP CAN BE ADDED LATER
###################
                        $i2=$fpr[$p];
                        $i22=$fpr[$p2];
                        $fine=0;
                        $longest = $bond_lst_l[$x][$i2];
                        if ($bond_lst_l[$x][$i2] < $bond_lst_l[$x][$i22] ) {$longest = $bond_lst_l[$x][$i22]};
                        if ($pr_mat[$p2][$x] == $pr_mat[$p][$x] && $longest > 1)
                        {
                            #printf "Current mess.\n";
                            # MAKE 2 DUMMY LISTS - DL1, DL2, INITIALISE LISTS TO SAME LENGTH
                            $dll = $bond_lst_l[$x][$i22];
                            if ( $bond_lst_l[$x][$i2] > $bond_lst_l[$x][$i22] )
                            {
                                $dll = $bond_lst_l[$x][$i2];
                            }
                            for ($dj = 1; $dj <= $nhat ; $dj++)     # sort will reorder everythign in this array, so it must be 0 for all possible length - alternatively if the sort function is truncated it would also work
                            {
                                $dl1[$dj]=0;
                                $dl2[$dj]=0;
                            }
   #                        FILL DUMMY LIST
                            for ($j=1 ; $j <= $bond_lst_l[$x][$i2] ; $j++)       # how many entries in x-bond list - go through each
                            {
                                $i3=$bond_list_ptr[$x][$i2][$j];
                                #printf"Blind luck  %i %i %i %i\n",   $lmat[1][$i3], $lmat[10][$i3], $longest, $dll;
                                $dl1[$j]=$lmat[10][$i3];
                            }
                            for ($j2=1 ; $j2 <= $bond_lst_l[$x][$i22] ; $j2++)       # how many entries in x-bond list - go through each
                            {
                                $i4=$bond_list_ptr[$x][$i22][$j2];
                                 #printf"Blind2 luck %i %i\n",   $lmat[1][$i4], $lmat[10][$i4];
                                $dl2[$j2]=$lmat[10][$i4];
                            }
#                           SORT DUMMY LISTS
                            @dl1s = sort { $b <=> $a } @dl1;
                            @dl2s = sort { $b <=> $a } @dl2;
#                           COMPARE DUMMY LISTS - IF THERE IS A WINNER RECORD THIS AND NOTE THAT THERE WAS A WINNER IN VARIABLE $FINE=1 AND EXIT
                            $fine=0;
                            for ($dj = 1; $dj <= $dll ; $dj++)
                            {
                               #printf "Comparing %i %i.\n", $dl1s[$dj] , $dl2s[$dj];

                                if ($dl1s[$dj] > $dl2s[$dj])
                                {
                                    #printf "i2 %i %i wins over %i. \n", $lmat[1][$i2], $dl1s[$dj] , $dl2s[$dj];
                                    $fine=1;
                                    last;
                                    #$dl1 -> $i2 is $p,
                                    #$dl2 -> $i22 is $p2
                                    
                                }
                                if ($dl1s[$dj] < $dl2s[$dj])
                                {
                                    #printf "i22 %i %i wins over %i.\n", $lmat[1][$i22], , $dl2s[$dj] , $dl1s[$dj];
                                    $temp_ch            = $cord[$p];
                                    $cord[$p]           = $cord[$p2];
                                    $cord[$p2]          = $temp_ch;

                                    for ($x2 = 0; $x2 < $nhat ; $x2++)
                                    {
                                        $temp_val           = $pr_mat[$p][$x2];
                                        $pr_mat[$p][$x2]    = $pr_mat[$p2][$x2];
                                        $pr_mat[$p2][$x2]   = $temp_val;
                                    }
                                    $fine=1;
                                    last;
                                }
                            }
                            if ($fine==1) {last} ;  # no need to keep looking at this bond for current p and p2 as it is resolved above
                            ############### start CIP R/S #################
                            if ($fine==0 )
                            {
                                #printf "not resolved %i %i %i %i %i\n",  $lmat[1][$fpr[$p]],  $lmat[1][$fpr[$p2]], $lmat[1][$i], $lmat[9][$i], $bond_lst_l[1][$i] ;
                                $found_4=0;

                                if ($bond_lst_l[1][$i] == 4)  {$found_4=1} ;    # four heteroatoms

                                if ($bond_lst_l[1][$i] >= 3)                    # three heteroatoms
                                {
                                    for ($xi=1 ; $xi <=$numlin ; $xi++)         #   print hydrogens
                                    {
                                        @fields = split (' ', $fil_c[$xi]);
                                        if ($xi > $atstart)
                                        {
                                            if($fields[0] == $lmat[1][$i])
                                            {
                                                @at4 = ($fields[2], $fields[3], $fields[4]);
                                                if($fields[5] =~ 'C.3')
                                                {
                                                    $found_4=1;
                                                }
                                            }
                                            if($fields[0] == $lmat[9][$i] && $lmat[9][$i] > 0)
                                            {
                                                #print $lmat[9][$i] $fields0;
                                                        @at3 = ($fields[2], $fields[3], $fields[4]);
                                            }
                                            if($fields[0] == $at_n && $lmat[9][$i] == 0)
                                            {
                                                #print $lmat[9][$i] $fields0;
                                                        @at3 = ($fields[2], $fields[3], $fields[4]);
                                            }

                                            if($fields[0] == $lmat[1][$fpr[$p]])
                                            {
                                                        @at1 = ($fields[2], $fields[3], $fields[4]);
                                            }
                                            if($fields[0] == $lmat[1][$fpr[$p2]])
                                            {
                                                        @at2 = ($fields[2], $fields[3], $fields[4]);
                                            }
                                            
                                        }
                                    }
                                }

                                    ##################
                                 if ($found_4==1)   # 3 bound heteroatoms - > sp2 or sp3, if no hydrogen then sp2
                                 {
                                     &getrs(@at1, @at2, @at3, @at4, $det);
                                     #printf "%f\n", $det;
                                     #printf "%f\n", $det;
                                     if ($det < 1)
                                     {
                                         #printf "Clockwise. order is: %i %i. \n", $lmat[1][$fpr[$p]],  $lmat[1][$fpr[$p2]];
                                         last;
                                     }
                                     if ($det > 1)
                                     {
                                         #printf "Anti-clockwise. order is: %i %i. \n", $lmat[1][$fpr[$p2]],  $lmat[1][$fpr[$p]];
                                         $temp_ch            = $cord[$p];
                                         $cord[$p]           = $cord[$p2];
                                         $cord[$p2]          = $temp_ch;

                                         for ($x2 = 0; $x2 < $nhat ; $x2++)
                                         {
                                             $temp_val           = $pr_mat[$p][$x2];
                                             $pr_mat[$p][$x2]    = $pr_mat[$p2][$x2];
                                             $pr_mat[$p2][$x2]   = $temp_val;
                                         }
                                         last;
                                     }
                                     #printf "%f %f %f\n", @at1;
                                     #printf "%f %f %f\n", @at2;
                                     #printf "%f %f %f\n", @at3;
                                     #printf "%f %f %f\n", @at4;
                                 }
                            }
                                ################### END CIP RULE
                            
                        }
                        if ($pr_mat[$p2][$x] > $pr_mat[$p][$x])
                        {
                            $temp_ch            = $cord[$p];
                            $cord[$p]           = $cord[$p2];
                            $cord[$p2]          = $temp_ch;

                            for ($x2 = 0; $x2 < $nhat ; $x2++)
                            {
                                $temp_val           = $pr_mat[$p][$x2];
                                $pr_mat[$p][$x2]    = $pr_mat[$p2][$x2];
                                $pr_mat[$p2][$x2]   = $temp_val;
                            }
                            last;
                        }
                        if ($pr_mat[$p2][$x] < $pr_mat[$p][$x])
                        {
                            last;
                        }
                    }


                }
            }
            for ($p = 1 ; $p <= $addchain ; $p++)       # PRINTING ONLY
            {
                #printf "Chain order is: %i %i\n", $p, $cord[$p];
                #printf "First atom in chain %i is: %i with priority: %i - chain ID: %i\n", $p, $lmat[1][$fpr[$p]], $cord[$p], $c_chains+($addchain-$cord[$p]);
                    #$c_chains=$c_chains+($addchain-$cord[$p]);
            }
############## ABOVE MATRIX CREATED
#############       END ASSIGNING PRIORITY  ##########
    
        }
    
    
#    $lmat[8][$i]=$t_chains;      #ca chain will be chain 1
    for ($j=$i+1 ; $j < $nhat ; $j++)
    {
        
        if ($lmat[1][$i] == $lmat[9][$j])   # This atom is connected to $i
        {

            $lmat[8][$j]=$lmat[8][$i];      #same chain id as what it is connected to.
            $lmat[11][$j]=$lmat[11][$i]+1;
            if ($addchain >= 1)
            {

                $c_chains=$addchain-1;
                for ($p = 1 ; $p <= $addchain ; $p++)
                {
                    if ($fpr[$p] == $j)
                    {
                        if ($cord[$p] == 1)
                        {
                            $lmat[8][$j]=$lmat[8][$i];
                        }else{
                            $lmat[8][$j] = $t_chains+$cord[$p]-1;
                        }
                            
                    #printf "Branch: %i %i %i %i %i\n", $lmat[1][$i], $lmat[1][$j], $lmat[8][$j], $t_chains, $cord[$p];
                    }
                }
            }
        }
        
    }
}


###########     create chains       #########
for ($j=1 ; $j <= $t_chains ; $j++)
{
    $cxl[$j]=0;
    $cxp[$j]=0;
    $cxptr[$j]=0;

    for ($i=1 ; $i < $nhat ; $i++)
    {
        if ($lmat[8][$i] == $j)
        {
                $k=$cxl[$j];
                $cx[$j][$k]=$lmat[1][$i];
                $cx_1[$j][$k]=$i;
                $cxl[$j]++;
        }
    }
#    for ($k = 0 ; $k < $cxl[$j] ; $k++)
#    {
#            printf "%i %i ", $cx[$j][$k], $cx_1[$j][$k];
#    }
#    printf "\n";
}
for ($j=1 ; $j <= $t_chains ; $j++)         ######### for printing only
{
    #printf "CHAIN %i: ", $j;
    for ($k = 0 ; $k < $cxl[$j] ; $k++)
    {
    #        printf "%i ", $cx[$j][$k]; #, $cx_1[$j][$k];
    }
    #printf "\n";
}

#   order based on length
for ($i = 1 ; $i < $chains ; $i++)
{
    for ($j = $i+1 ; $j <= $chains ; $j++)
    {
        if ($cxl[$j] > $cxl[$i])
        {
            $tvl1       = $cxl[$i];
            $tprt       = $cxp[$i];
            for ($k = 0 ; $k < $cxl[$i] ; $k++)
            {
                $tvcx1[$k]  = $cx[$i][$k];
                $tvcx5[$k]  = $cx_1[$i][$k];
            }
            for ($k = 0 ; $k < $cxl[$j] ; $k++)
            {
                   $cx[$i][$k]     =    $cx[$j][$k];
                $cxptr[$i][$k]     =  $cx_1[$j][$k];
            }
            $cxl[$i]        = $cxl[$j];
            $cxp[$i]        = $cxp[$j];

            $cxl[$j]        = $tvl1;
            $cxp[$j]        = $tprt;
            for ($k = 0 ; $k < $tvl1 ; $k++)
            {
                   $cx[$j][$k]     = $tvcx1[$k];
                 $cx_1[$j][$k]     = $tvcx5[$k];
            }

        }
        
    }
}
#$cxp[2]=5;

######################################
#   It's a bit unclear if at a branch point you simply continue to name based on priority or make new branches with new priorities or identify which is the main chain based on the priority of that chain - i think this is probably what it will be - and because it has the lowest priority
#   It's clear that it follows taht the naming follows the chain structure. i.e. tryptophan does not have a z1!
#   If same priority at atom +1 -> number of connected atoms at +1 -> priority at atom +2 -> if same number of connected atoms at +2
#
##########
####### HOW MANY POSITIONS ARE THERE? ############
$maxlmat=0;
for ($i=1 ; $i < $nhat ; $i++)
{
#    printf "Atom %i Position: %s Dist: %i\n", $lmat[1][$i], $gnam[$lmat[11][$i]], $lmat[11][$i];
    if ($lmat[11][$i] > $maxlmat)
    {
        $maxlmat = $lmat[11][$i];
    }
}
for ($x=0 ; $x < $maxlmat ; $x++)   # initialise padded chain vector
{
    for ($j=1 ; $j <= $t_chains ; $j++)
    {
            $p_cx[$j][$x]=0;
    }
    $ppos[$x]=0;
}
#printf "Longest: %i Position: %s \n", $maxlmat, $gnam[$maxlmat];
########### FILL PADDED VECTOR
for ($j=1 ; $j <= $t_chains ; $j++)
{
    for ($k = 0 ; $k < $cxl[$j] ; $k++)
    {
            $t1=$cx_1[$j][$k];
            $t2=$lmat[11][$t1];
#            printf "%i %i ", $cx[$j][$k], $t2;
            $p_cx[$j][$t2]=$cx[$j][$k];
    }
    #printf "\n";
}
################ COMPARE POSITIONALLY
for ($x=0 ; $x <= $maxlmat ; $x++)
{
    $count_matl=0;
    for ($j=1 ; $j <= $t_chains ; $j++)
    {
#            printf "%i\t", $p_cx[$j][$x];
                if ($p_cx[$j][$x] > 0)
                {
                   $count_matl++;
                }
    }
    $c_ch_at[$x]=$count_matl;
#    printf "\n";

}
#exit;


for ($i=1 ; $i < $nhat ; $i++)
{
    $lmat[13][$i]='0';
#    $lmat[15][$k]=0;
    $lmat[16][$k]=0;
    #printf "%i.\n", length($lmat[2][$i]);
}
printf "Updated atom names.\n";

for ($i=1 ; $i < $nhat ; $i++)
{
 # if this is main chain and only entry
    if ($lmat[8][$i] == 1 && $c_ch_at[$lmat[11][$i]] == 1)
    {
        $lmat[13][$i] = sprintf("%s%s ", $lmat[2][$i], $gnam[$lmat[11][$i]]);
    }else{
        $lmat[13][$i] = sprintf("%s%s%i", $lmat[2][$i], $gnam[$lmat[11][$i]], $lmat[8][$i]);
    }
#    if ()
    printf "%i %s \n", $lmat[1][$i], $lmat[13][$i];
}
#exit;


##############   CODE FOR PRINTING HYDROGENS   #################
$l=0;
$allnh=0;
for ($i=1 ; $i <=$numlin ; $i++)        #   print hydrogens
{
    @fields = split (' ', $fil_c[$i]);
    if ($i > $atstart)
    {
        $hlen=length($fields[5]);
        if($fields[5] =~ 'H' && $hlen == 1)
        {
            #printf "%s \n", $fil_c[$i];

            $l++;
            $hmat[1][$l]=$fields[0];
            $hmat[4][$l]=0;
            #$hmat[2][$l]=$fields[8];
            #$hmat[7][$l]=$fields[12];       #   pseudo atoms
            $hmat[8][$l]=0;
            $hmat[21][$l]=0;
                    for ($j = 1 ; $j <= $numbnd ; $j++)
                    {
                        #printf "%s %s %s\n", $temp_b1[$j], $temp_b2[$j], $hmat[1][$l];

                        if ($temp_b1[$j] == $hmat[1][$l])
                        {
                         #   printf "xxx%i %i %i\n", $temp_b1[$j], $temp_b2[$j], $at_ca;
                            $hmat[2][$l]=$temp_b2[$j];
                            last;
                        }
                        if ($temp_b2[$j] == $hmat[1][$l])
                        {
                          #  printf "ytyy%i %i %i\n", $temp_b1[$j], $temp_b2[$j], $at_ca;
                            $hmat[2][$l]=$temp_b1[$j];
                            last;
                        }
                    }

            for ($k=1 ; $k < $nhat ; $k++)
             {
                if($lmat[1][$k] == $hmat[2][$l])
                {
                    $hmat[3][$l]=$lmat[11][$k];         # so now I know it is hb or hg etc.
                    $lmat[16][$k]++;                    #lmat 16 = number of H's attached to that carbon total
                    $hmat[5][$l]=$lmat[16][$k];
                    #printf "%i %i\n", $lmat[1][$k], $hmat[5][$l];

                }
            }
        if ($fields[8] == 1 && $nterm == 1)
        {
         $hmat[21][$l]++;
         $allnh=$hmat[21][$l];
        }
        }
    }
}

$hhat=$l;

for ($i=1 ; $i <=$hhat ; $i++)
{
    for ($k=1 ; $k < $nhat ; $k++)
    {
        if ($lmat[1][$k] == $hmat[2][$i] && $hmat[4][$i] == 0)
        {
            $hmat[4][$i]=$lmat[16][$k];          #lmat 16 = number of H's attached to that carbon total - also how many h's there are
            #printf "3at no: %i carb: %i \n", $lmat[1][$k], $hmat[2][$i];
            $hmat[20][$i]=$lmat[8][$k];          #same chain, used later for printing
            if ($lmat[8][$k] == 1 && $c_ch_at[$lmat[11][$k]] == 1)
            {
                    if ($hmat[4][$i] == 1)
                    {
                        $hmat[6][$i] = sprintf("H%s  ", $gnam[$lmat[11][$k]]);
                    }elsif ($hmat[4][$i] == 2 && $lmat[10][$k] == 6)        # for carbons it is HX2 and HX3 for some reason.
                    {
                        ####### CIP FOR HYDROGENS HX2 AND HX3 #############
                        $anti = 0;
                        
                        if ($hmat[5][$i] < $hmat[4][$i])    # if this is the first of a pair
                        {
                            $lmat[30][$k]=$hmat[1][$i];
                        }
                        if ($hmat[5][$i] == $hmat[4][$i])    # if this is the first of a pair
                        {
                            #                            $lmat[30][$k]=$hmat[1][$i];
                            #printf "central: %i back: %i  h1: %i h2:%i\n", $lmat[1][$k], $lmat[9][$k], $hmat[1][$i], $lmat[30][$k];
                                    
                            for ($xi=1 ; $xi <=$numlin ; $xi++)         #   print hydrogens
                            {
                                @fields = split (' ', $fil_c[$xi]);
                                if ($xi > $atstart)
                                {
                                    #printf "%i %f %f\n",  $fields[0], $fields[6], $fields[7];
                                    
                                    if($fields[0] == $lmat[1][$k])        # PUT SUBSTIUENT WHERE THE CENTRAL ATOM IS AND GIVE LOWEST PRIORITY
                                    {
                                        @at4 = ($fields[5], $fields[6], $fields[7]);
                                    }
                                    if($fields[0] == $at_n && $lmat[9][$i] == 0)    # first nitrogen you hit - shoudl be 1.
                                    {
                                        @at3 = ($fields[5], $fields[6], $fields[7]);
                                    }
                                    
                                    if($fields[0] == $lmat[9][$i] && $lmat[9][$i] > 0)
                                    {
                                        #print $lmat[9][$i] $fields0;
                                        @at3 = ($fields[5], $fields[6], $fields[7]);
                                    }
                                    if($fields[0] == $lmat[30][$k])        # first hydrogen atom encountered on that carbon
                                    {
                                        @at2 = ($fields[5], $fields[6], $fields[7]);
                                    }
                                    if($fields[0] == $hmat[1][$i])
                                    {
                                        @at1 = ($fields[5], $fields[6], $fields[7]);
                                    }
                                }
                            }
                            
                            &getrs(@at1, @at2, @at3, @at4, $det);
                            if ($det < 1)
                            {
                                #printf "Clockwise.\n";
                            }
                            if ($det > 1)       # should only be here on second atom of 2
                            {
                                #printf "Anti-clockwise. \n";
                                $anti=1;
                            }
                            
                            
                        }
                        ##################    END   for CIP   ##################

                        $hmat[6][$i] = sprintf("H%s%i ", $gnam[$lmat[11][$k]], $hmat[5][$i]+1);
                        ##################     for CIP2   ###################
                        if ($anti == 0 && $hmat[5][$i] < $hmat[4][$i])
                        {
                            $lmat[31][$k]=$hmat[6][$i];
                            $lmat[32][$k]=$i;
                        }
                        if ($anti == 1 && $hmat[5][$i] == $hmat[4][$i])
                        {
                            $hmat[6][$lmat[32][$k]]=$hmat[6][$i];
                            $hmat[6][$i]=$lmat[31][$k];
                            #printf "swap\n";
                        }
                        ##################    END   for CIP2   ##################
                    }else{
                        $hmat[6][$i] = sprintf("H%s%i ", $gnam[$lmat[11][$k]], $hmat[5][$i]);
                    }
            }else{
                #$lmat[13][$i] = s
                if ($hmat[4][$i] == 1)
                {
                    $hmat[6][$i] = sprintf("H%s%i ", $gnam[$lmat[11][$k]], $lmat[8][$k]);
                }elsif ($hmat[4][$i] == 2  && $lmat[10][$k] == 6)       # for carbons it is HX2 and HX3 for some reason.
                {
                    $anti = 0;
                     ####### CIP FOR HYDROGENS HX2 AND HX3 #############
                     if ($hmat[5][$i] < $hmat[4][$i])    # if this is the first of a pair
                     {
                         $lmat[30][$k]=$hmat[1][$i];
                     }
                     if ($hmat[5][$i] == $hmat[4][$i])    # if this is the first of a pair
                     {
                         #                            $lmat[30][$k]=$hmat[1][$i];
                         #printf "%i %i  %i %i\n", $lmat[1][$k], $lmat[9][$k], $hmat[1][$i], $lmat[30][$k];
                         for ($xi=1 ; $xi <=$numlin ; $xi++)         #   print hydrogens
                         {
                             @fields = split (' ', $fil_c[$xi]);
                             if ($xi > $atstart)
                             {
                                 #printf "%i %f %f\n",  $fields[0], $fields[6], $fields[7];
                                 
                                 if($fields[0] == $lmat[1][$k])        # PUT SUBSTIUENT WHERE THE CENTRAL ATOM IS AND GIVE LOWEST PRIORITY
                                 {
                                     @at4 = ($fields[5], $fields[6], $fields[7]);
                                 }
                                 if($fields[0] == $at_n && $lmat[9][$i] == 0)    # first nitrogen you hit - shoudl be 1.
                                 {
                                     @at3 = ($fields[5], $fields[6], $fields[7]);
                                 }
                                 
                                 if($fields[0] == $lmat[9][$i] && $lmat[9][$i] > 0)
                                 {
                                     #print $lmat[9][$i] $fields0;
                                     @at3 = ($fields[5], $fields[6], $fields[7]);
                                 }
                                 if($fields[0] == $lmat[30][$k])        # first hydrogen atom encountered on that carbon
                                 {
                                     @at2 = ($fields[5], $fields[6], $fields[7]);
                                 }
                                 if($fields[0] == $hmat[1][$i])
                                 {
                                     @at1 = ($fields[5], $fields[6], $fields[7]);
                                 }
                             }
                         }
                         
                         &getrs(@at1, @at2, @at3, @at4, $det);
                         if ($det < 1)
                         {
                             #printf "Clockwise.\n";
                         }
                         if ($det > 1)       # should only be here on second atom of 2
                         {
                             #printf "Anti-clockwise. \n";
                             $anti=1;
                         }
                         
                         
                     }
                     ##################    END   for CIP   ##################
                    $hmat[6][$i] = sprintf("H%s%i%i", $gnam[$lmat[11][$k]], $lmat[8][$k], $hmat[5][$i]+1);
                    ##################     for CIP2   ###################
                    if ($anti == 0 && $hmat[5][$i] < $hmat[4][$i])
                    {
                        $lmat[31][$k]=$hmat[6][$i];
                        $lmat[32][$k]=$i;
                    }
                    if ($anti == 1 && $hmat[5][$i] == $hmat[4][$i])
                    {
                        $hmat[6][$lmat[32][$k]]=$hmat[6][$i];
                        $hmat[6][$i]=$lmat[31][$k];
                        #printf "swap\n";
                    }
                    ##################    END   for CIP2   ##################
                }else{
                    $hmat[6][$i] = sprintf("H%s%i%i", $gnam[$lmat[11][$k]], $lmat[8][$k], $hmat[5][$i]);
                }
            }
        }
    }
    if ($hmat[2][$i] == 3 && $nterm == 0)
    {
        $hmat[6][$i] = sprintf("H   ");     #special case for amide H
        $hmat[20][$i]=1;
    }
    if ($hmat[2][$i] == 1 && $nterm == 1)
    {
        if ($allnh == 1)
        {
        $hmat[6][$i] = sprintf("H   ");     #special case for amide H
        }else{
        $hmat[6][$i] = sprintf("H%i  "), $hmat[21][$i];     #special case for amide H
        }
        $hmat[20][$i]=1;
    }
}
for ($i=1 ; $i <=$hhat ; $i++)        # get pseudoatoms
{
    if ($hmat[5][$i] > 0)
    {
            printf "%i %s\n", $hmat[1][$i], $hmat[6][$i];
    }
}

############### PRINT #############



#$mol2 = sprintf("block_%s", $inputfile);

open(OUT, '>', $outputfile) || die "Can't open output file";
#        printf OUT "%s",$head[$i];

$endb=0;
for ($j = 1 ; $j < $numlin ;  $j++)
{
    @fields = split (' ', $fil_c[$j]);

     if ($fil2_c[$j] == 3)
     {
            $endb=0;
            for ($i=1 ; $i <= $hhat ; $i++)      ### just for printing
            {
                if ($hmat[1][$i] == $fields[0] && $hmat[5][$i] > 0)
                {
                    $hmat[6][$i]=~ s/^\s+|\s+$//g;
                    #printf "%i\n", length($lmat[13][$i]);
                    $hmat[6][$i]=sprintf("%4s", $hmat[6][$i]);
                    substr($fil_c[$j],8,4)=$hmat[6][$i];
                    printf OUT "%s\n", $fil_c[$j];
#                printf "%i %s\n",$hmat[1][$i], $hmat[6][$i];
                $endb=1;
                last;
                }
            }
            for ($i=1 ; $i < $nhat ; $i++)      ### just for printing
            {
             if ($lmat[1][$i] == $fields[0])
             {
                 #$lmat[13][$i]=sprintf("BC");
                 #printf "%i\n", length($lmat[13][$i]);
                 $lmat[13][$i]=~ s/^\s+|\s+$//g;
                 #printf "%i\n", length($lmat[13][$i]);
                 $lmat[13][$i]=sprintf("%4s", $lmat[13][$i]);
                 substr($fil_c[$j],8,4)=$lmat[13][$i];
                 printf OUT "%s\n", $fil_c[$j];
             $endb=1;
             last;
             }
            }
         if ($endb==0){printf OUT "%s\n", $fil_c[$j]}
     }else{
     printf OUT "%s\n", $fil_c[$j];
     }
}
close OUT;







############### PRINT #############



exit;




##############   CODE FOR PRINTING HYDROGENS END  #################




sub getatno($attyp,$atno)
{
## NEED REVERSE OF THIS
    $atno=0;
    if ($attyp =~ 'C') {$atno = 6};
    if ($attyp =~ 'N') {$atno = 7};
    if ($attyp =~ 'O') {$atno = 8};
    if ($attyp =~ 'F') {$atno = 9};
    if ($attyp =~ 'S') {$atno = 16};
    if ($attyp =~ 'Si') {$atno = 14};
    if ($attyp =~ 'P') {$atno = 15};
    if ($attyp =~ 'Cl') {$atno = 17};
    if ($attyp =~ 'Br') {$atno = 35};
    if ($attyp =~ 'I') {$atno = 53};
    if ($atno == 0)
    {
        $atno == 99;
         printf "%s not defined. Atom will be given atom number %i.\n", $attyp, $atno;
    }
   # printf "%s %i\n", $attyp, $atno;

## to convert between 3 and 1 letter codes
}

sub getattyp($atno,$attyp)
{
## NEED REVERSE OF THIS
    $attyp = ' ';

    if ( $atno == 6 ) { $attyp = 'C' } ;
    if ( $atno == 7 ) { $attyp = 'N' } ;
    if ( $atno == 8 ) { $attyp = 'O' } ;
    if ( $atno == 9 ) { $attyp = 'F' } ;
    if ( $atno == 16 ) { $attyp = 'S' } ;
    if ( $atno == 14 ) { $attyp = 'Si' } ;
    if ( $atno == 15 ) { $attyp = 'P' } ;
    if ( $atno == 17 ) { $attyp = 'Cl' } ;
    if ( $atno == 35 ) { $attyp = 'Br' } ;
    if ( $atno == 53 ) { $attyp = 'I' } ;

}
    
    sub getrs(@at1, @at2, @at3, @at4)
    {


    $A11=$at1[0];
    $A12=$at1[1];
    $A13=$at1[2];
    $A14=1;
    $A21=$at2[0];
    $A22=$at2[1];
    $A23=$at2[2];
    $A24=1;
    $A31=$at3[0];
    $A32=$at3[1];
    $A33=$at3[2];
    $A34=1;
    $A41=$at4[0];
    $A42=$at4[1];
    $A43=$at4[2];
    $A44=1;
    # determinant
    $det = $A11*$A22*$A33*$A44 + $A11*$A23*$A34*$A42 + $A11*$A24*$A32*$A43
          -$A11*$A24*$A33*$A42 - $A11*$A23*$A32*$A44 - $A11*$A22*$A34*$A43
          -$A12*$A21*$A33*$A44 - $A13*$A21*$A34*$A42 - $A14*$A21*$A32*$A43
          +$A14*$A21*$A33*$A42 + $A13*$A21*$A32*$A44 + $A12*$A21*$A34*$A43
          +$A12*$A23*$A31*$A44 + $A13*$A24*$A31*$A42 + $A14*$A22*$A31*$A43
          -$A14*$A23*$A31*$A42 - $A13*$A22*$A31*$A44 - $A12*$A24*$A31*$A43
          -$A12*$A23*$A34*$A41 - $A13*$A24*$A32*$A41 - $A14*$A22*$A33*$A41
          +$A14*$A23*$A32*$A41 + $A13*$A22*$A34*$A41 + $A12*$A24*$A33*$A41;
    }


exit;





