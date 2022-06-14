#!/usr/bin/perl
#
#Program to rename MOL2 file according to IUPAC
#
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);

#READ OPTIONS ****************************
GetOptions('ca=s' => \$ca,
           'n=s' => \$n,
           'c=s' => \$c,
           'i=s' => \$inputfile,
           'o=s' => \$outputfile
			);
# -c C6 -x1 C7 -x2 N2
#READ LIB FILE **************************
open(LIB, $inputfile) || die "Error: Could not open input library file.\nUsage\n ./fixnom_mini.pl -i input.mol2 -o output.mol2 -ca <CA atom name> -c <C atom name> -n <N atom name>\n ";		# Open the input file
@rlib = <LIB>;										# Read it into an array
close(LIB);
# ***********************
# READ GREEK NAMES
my @gnam  =('A', 'B', 'G', 'D', 'E', 'Z', 'H', 'Q', 'I', 'K', 'L', 'M', 'N', 'X', 'O', 'P', 'R', 'S', 'T', 'Y', 'F', 'C', 'U', 'W');
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

            #printf "%s %s %s %s\n", $fields[0], $fields[1], $fields[5], $temp_c2[$l] ;
            # printf "%s %s %s %s %s\n", $c, $n, $ca, $fields[1], $temp_c2[$l] ;

            if ($ca =~ $temp_c2[$l]  &&  length($ca) == length($temp_c2[$l]))
            {
                $at_ca=$temp_c1[$l];
             #   printf "ca: %s\n", $at_ca;
            }
            if ($c =~ $temp_c2[$l]  &&  length($c) == length($temp_c2[$l]))
            {
                $at_c=$temp_c1[$l];
             #   printf "c: %s\n", $at_c;
            }
            if ($n =~ $temp_c2[$l]  &&  length($n) == length($temp_c2[$l]))
            {
                $at_n=$temp_c1[$l];
             #   printf "n: %s\n", $at_n;
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
# temp_c1 goes from 1 to numhat and has the atom ID.
# temp_b1 goes from 1 to numbnd and has the bonds beteween atoms.
#################################get NH AND CO ########################
$at_hn=$at_co=$at_tco=$at_hn2=$at_hn3=-1;
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
                    if ($at_co > 0) # CO already found so this second O must be part of C-terminus.
                    {
                        $at_tco=$temp_b2[$j];
                    }else{
                    $at_co = $temp_b2[$j] ;
                    }
                }
                if($temp_b2[$j] == $temp_c1[$i] && $temp_c3[$temp_b1[$j]] =~ 'O')
                {
                    if ($at_co > 0) # CO already found so this second O must be part of C-terminus.
                    {
                        $at_tco=$temp_b1[$j];
                    }else{
                    $at_co = $temp_b1[$j] ;
                    }
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
                    if ($at_hn > 0 && $at_hn2 > 0) # NH already found so this second/thrid NH  must be part of N-terminus. only goes to NH3
                    {
                        $at_hn3=$temp_b2[$j];
                    }elsif($at_hn > 0 && $at_hn2 < 0)
                    {
                        $at_hn2=$temp_b2[$j];
                    }elsif($at_hn < 0 )
                    {
                        $at_hn = $temp_b2[$j] ;
                    }
                }
                if($temp_b2[$j] == $temp_c1[$i] && $temp_c3[$temp_b1[$j]] =~ 'H')
                {
                    if ($at_hn > 0 && $at_hn2 > 0) # NH already found so this second/thrid NH  must be part of N-terminus. only goes to NH3
                    {
                        $at_hn3=$temp_b1[$j];
                    }elsif($at_hn > 0 && $at_hn2 < 0)
                    {
                        $at_hn2=$temp_b1[$j];
                    }elsif($at_hn < 0 )
                    {
                        $at_hn = $temp_b1[$j] ;
                    }
#                    $at_hn = $temp_b1[$j] ;
                }
            }
        }
    }
}

# above code identifies the HN hydrogen and teh CO oxygen. - at_nh and at_co

################################# CUT OUT CO SIDE OF THINGS ######################
# NEED TO reduce matrix further to cut out N/C connections.
# assumes index goes form 1 to end.
# Sort so that CO is on top
for ($i=1 ; $i <= $numhat ; $i++)
{
    if($temp_c1[$i] == $at_c && $temp_c1[1] != $temp_c1[$i])
    {
        ($temp_c1[1], $temp_c1[$i]) = ($temp_c1[$i], $temp_c1[1]);
        ($temp_c2[1], $temp_c2[$i]) = ($temp_c2[$i], $temp_c2[1]);
        ($temp_c3[1], $temp_c3[$i]) = ($temp_c3[$i], $temp_c3[1]);
        ($temp_c4[1], $temp_c4[$i]) = ($temp_c4[$i], $temp_c4[1]);

    }
}

#printf "oxygen %i %i %i\n", $at_c, $at_n, $at_ca;
$temp_c4[1]=1;      # above C' is set to 1 so this is C'
$t1=1;
for ($i=1 ; $i <= $numhat ; $i++)
{
    if($temp_c4[$i] == 1)           # only true when looking at C'
    {
        for ($j = 1 ; $j <= $numbnd ; $j++)
        {
            #printf "A: %i %i %i\n", $temp_c1[$i], $temp_b1[$j], $at_ca;

            unless ($temp_b1[$j] == $at_ca || $temp_b2[$j] == $at_ca)
            {
                if($temp_b1[$j] == $temp_c1[$i] )
                {
                    #printf "xxx%i %s %i\n", $temp_b1[$j], $temp_c1[$i], $i;
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
    for ($j=$i+1 ; $j < $numhat ; $j++)
    {
        for ($k=$j+1 ; $k <= $numhat ; $k++)
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
#for ($i=1 ; $i <= $numhat ; $i++)
#{
#    print "$temp_c1[$i] $temp_c2[$i] $temp_c3[$i] $temp_c4[$i]\n";
#}
################################# CUT OUT N SIDE OF THINGS ######################
# NEED TO reduce matrix further to cut out N/C connections.
# assumes index goes form 1 to end.
# Sort so that N is on top
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

            unless ($temp_b1[$j] == $at_ca || $temp_b2[$j] == $at_ca)   # fails for proline
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
    for ($j=$i+1 ; $j < $numhat ; $j++)      ### just for printing
    {
        for ($k=$j+1 ; $k <= $numhat ; $k++)      ### just for printing
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

for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c1[$i] == $at_ca && $temp_c1[1] != $temp_c1[$i])
    {
        ($temp_c1[1], $temp_c1[$i]) = ($temp_c1[$i], $temp_c1[1]);
        ($temp_c2[1], $temp_c2[$i]) = ($temp_c2[$i], $temp_c2[1]);
        ($temp_c3[1], $temp_c3[$i]) = ($temp_c3[$i], $temp_c3[1]);
        ($temp_c4[1], $temp_c4[$i]) = ($temp_c4[$i], $temp_c4[1]);

    }
}

#printf "oxygen %i\n", $at_co;
$temp_c4[1]=-2;

for ($i=1 ; $i <= $numhat ; $i++)      ### just for printing
{
    if($temp_c4[$i] == -2)
    {
        for ($j = 1 ; $j <= $numbnd ; $j++)
        {
#            printf "%i %i %i\n", $temp_c1[$i], $temp_b1[$j], $at_ca;

#            unless ($temp_b1[$j] == $at_c || $temp_b2[$j] == $at_c || $temp_b2[$j] == $at_n || $temp_b2[$j] == $at_n)   # fails for proline
#            {
                if($temp_b1[$j] == $temp_c1[$i] )
                {
                    #printf "%i %s %i\n", $temp_b2[$j], $temp_c2[$temp_b2[$j]], $i;
                    #need to match b2 to atom list
                    for ($i2=$i+1 ; $i2 <= $numhat ; $i2++)      ### just for printing
                    {
                        if ($temp_c1[$i2] == $temp_b2[$j])
                        {
                            unless ($temp_c1[$i2] == $at_c || $temp_c1[$i2] == $at_n)
                            {
                             $temp_c4[$i2]=-2;
                            }
                        }
                    }
                        
                }
                if($temp_b2[$j] == $temp_c1[$i] )
                {

                    for ($i2=$i+1 ; $i2 <= $numhat ; $i2++)      ### just for printing
                    {
                        if ($temp_c1[$i2] == $temp_b1[$j])
                        {
                            unless ($temp_c1[$i2] == $at_c || $temp_c1[$i2] == $at_n)
                            {
                             $temp_c4[$i2]=-2;
                            }
                        }
                    }
                }
                
 #           }
        }
    }


    ## sort the array but record the ordered indexes.
    for ($j=$i+1 ; $j < $numhat ; $j++)      ### just for printing
    {
        for ($k=$j+1 ; $k <= $numhat ; $k++)      ### just for printing
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

for ($i=1 ; $i <= $numhat ; $i++)
{
    if ($temp_c4[$i] == -2) {$temp_c4[$i]=0};
#    print "E: $temp_c1[$i] $temp_c2[$i] $temp_c3[$i] $temp_c4[$i]\n";
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

my($filename, $directories, $suffix) = fileparse($inputfile);
$mol2 = sprintf("%s/block_%s", $directories, $filename);

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
#    printf "mid %i %i %i %i\n", $at_ca, $at_c, $nhat,  $lmat[1][10];

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
$dum=0;
#
#
$test=1;
if ($test == 1)
{
    
    for ($i = 1 ; $i <= $numbnd ; $i++)
    {
       # printf "OUT %i %i %s %i %i\n", $temp_b1[$i], $temp_b2[$i], $temp_b3[$i], $temp_c6[$temp_b1[$i]], $temp_c6[$temp_b2[$i]] ;
        if ($temp_b3[$i] =~ 'ar')
        {
         #   printf "%i %i %s %i %i\n", $temp_b1[$i], $temp_b2[$i], $temp_b3[$i], $temp_c6[$temp_b1[$i]], $temp_c6[$temp_b2[$i]] ;

            $temp_b3[$i]=2;
            $temp_c6[$temp_b1[$i]]++;
            $temp_c6[$temp_b2[$i]]++;       # contineu from here
            

        }
        if ($temp_b3[$i] > 1)
        {
           # printf "%i %i\n", $temp_b1[$i], $temp_b3[$i];
            # find the two atoms and add in the dummy for each
              #  printf "I'm here with %i %i\n", $temp_b2[$i], $temp_b1[$i]  ;

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
              #  printf "INSIDE %i %i %s %i %i\n", $temp_b1[$i], $temp_b2[$i], $temp_b3[$i], $temp_c6[$temp_b1[$i]], $temp_c6[$temp_b2[$i]] ;

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
                        $lmat[12][$nhat+$dum]=$lmat[1][$id1];      #PAIRWISE CONNECTION (BACK TO CA)
                     # dummy atom has 1 1-bond partner
                     #      $bond_lst_l[1][$nhat+$dum]=1;
                     # dummy atoms 1 bond partner is
                     #      $bond_list[1][$nhat+$dum][1]=$lmat[1][$id1];
                     
                     # the atom connected to dummy now has 1 more 1-bond
                           $bond_lst_l[1][$id1]++;
                           $b_l=$bond_lst_l[1][$id1];
                     # the atom connected to dummy has dummy in that position ... ptr updated later
                           $bond_list[1][$id1][$b_l]=$lmat[1][$nhat+$dum];
                   #  printf "dummy ATOM %i IS A %s CONNECTED TO ATOM %i WITH THE %i BOND. ch id %i p %i\n", $nhat+$dum, $lmat[2][$nhat+$dum], $lmat[9][$nhat+$dum], $l, $lmat[3][$nhat+$dum], $lmat[10][$nhat+$dum];

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
                           $lmat[12][$nhat+$dum]=$lmat[1][$id2];
                     # dummy atom has 1 1-bond partner
                     #      $bond_lst_l[1][$nhat+$dum]=1;
                     # dummy atoms 1 bond partner is
                     #      $bond_list[1][$nhat+$dum][1]=$lmat[1][$id2];
                     
                     # the atom connected to dummy now has 1 more 1-bond
                           $bond_lst_l[1][$id2]++;
                           $b_l=$bond_lst_l[1][$id2];
                     # the atom connected to dummy has dummy in that position ... ptr updated later
                           $bond_list[1][$id2][$b_l]=$lmat[1][$nhat+$dum];
                   # printf "dummy2 ATOM %i IS A %s CONNECTED TO ATOM %i WITH THE %i BOND. ch id %i p %i\n", $nhat+$dum, $lmat[1][$nhat+$dum], $lmat[9][$nhat+$dum], $l, $lmat[3][$nhat+$dum], $lmat[10][$nhat+$dum];

                    $dum++;
                    }
                
                }
            }
            
        }

    }
    $nhat=$nhat+$dum;

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
#
#    for ($i2 = 1 ; $i2 < $nhat ; $i2++)
#    {
##        printf "%i %i %i %s %i\n", $dum, $lmat[1][$i2], $lmat[10][$i2], $lmat[2][$i2], $lmat[9][$i2];
##    exit;
#    }
#

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
            #printf "%i ", $bond_list[1][$i][$j];
            $master_l[$i][$j] = $bond_list[1][$i][$j];   # first entry is itself
        }
    
    #printf "\n";
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
   # printf "%i %i %i %s %i\n", $dum, $lmat[1][$i], $lmat[10][$i], $lmat[2][$i], $lmat[9][$i];

    $bond_list[0][$i][1]=$lmat[1][$i];
    $bond_list_ptr[0][$i][1]=$i;
    $bond_lst_l[0][$i]=1;
    for ($x = 1; $x < $nhat ; $x++)
    {
        for ($j=1 ; $j <= $bond_lst_l[$x][$i] ; $j++)       # how many entries in x-bond list - go through each
        {
            $dist_ca=-1;
            if ($bond_list[$x][$i][$j] == $at_ca)
            {
                $dist_ca[$i]=$x;
            }
            #printf "id = %i, i = %i, j= %i, x = %i, list:  %i\n", $lmat[1][$i], $i, $j, $x, $bond_list[$x][$i][$j];
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
#printf "Modified: %i\n", $lmat[1][10];
$nhat=$nhat-$dum;
for ($i=1 ; $i < $nhat-1 ; $i++)
{
    $dumi[$i]=$i;
}


for ($dmi=1 ; $dmi < $nhat-1 ; $dmi++)
{
    $i=$dumi[$dmi];
####### NEED TO MAKE SURE I SORT FUTURE ATOMS ONCE CHAINS ARE BEING ASSIGNED
    #printf "xxxx: %i %i %i\n", $lmat[1][$i], $dist_ca[$i], $i;
    $isort = 0;

    for ($ix=1 ; $ix < $nhat-1 ; $ix++)
    {
        if($dist_ca[$ix] == $dist_ca[$i]+1)
        {
           
            #printf "id: %i, chain: %i, seen: %i, dist: %i, connected to: %i ", $lmat[1][$ix], $lmat[8][$ix], $lmat[11][$ix], $dist_ca[$ix], $lmat[9][$ix];
            $tsrtid[$isort]=$lmat[1][$ix];
            $tsrtpr[$isort]=0;
            $tsrtpt[$isort]=$ix;
            $tsrtpr[$isort]=$lmat[8][$ix];
            $isort++;
        }
    }
    for ($ix=0; $ix < $isort-1 ; $ix++)
    {
        for ($jx=$ix+1; $jx < $isort ; $jx++)
        {
            if($tsrtpr[$ix] > $tsrtpr[$jx])
            {
#                 ($first, $second) = ($second, $first);
                ($tsrtpr[$ix], $tsrtpr[$jx]) =  ($tsrtpr[$jx], $tsrtpr[$ix]);
                ($tsrtid[$ix], $tsrtid[$jx]) =  ($tsrtid[$jx], $tsrtid[$ix]);
                $sw1=$tsrtpt[$ix];
                $sw2=$tsrtpt[$jx];
                (  $dumi[$sw1],   $dumi[$sw2]) =  (  $dumi[$sw2],   $dumi[$sw1]);
                #printf "swap atoms %i and %i, %i %i\n", $tsrtid[$ix], $tsrtid[$jx], $dumi[$sw1],   $dumi[$sw2];
            }
        }
    }
### END    ####### NEED TO MAKE SURE I SORT FUTURE ATOMS ONCE CHAINS ARE BEING ASSIGNED

    
#    $attyp = $lmat[2][$i];
#    getatno($attyp, $atno);
#    $lmat[10][$i]=$atno;
    $t_chains=$t_chains+$c_chains;
    $c_chains=0;
    $found=0;
    $addchain=0;
    $fpr_i[$i]=0;

############        ASSIGNING PRIORITY  ############
                if ($bond_lst_l[1][$i]>2 || ($bond_lst_l[1][$i] == 2 && $lmat[1][$i] == $at_ca) )               # if more than 2 attached heteroatmos, then this is a branch
                {
                    for ($j2=1 ; $j2 <= $bond_lst_l[1][$i] ; $j2++)         # check how many branching atoms > $i - > if less than $i then not branch
                    {
                        if ($bond_list_ptr[1][$i][$j2] > $i && $lmat[1][$bond_list_ptr[1][$i][$j2]] < 1000)
                        #if ($bond_list_ptr[1][$i][$j2] > $i )                # since things are sorted above if the id is less than current ID then it can be ignored
                        {
                            if ($dist_ca[$bond_list_ptr[1][$i][$j2]] > $dist_ca[$i])
                            {
                            $addchain++;                                    # if more than 1 you really do have a new branch - addchain tells you how many chains needed
                                # probably should not increment if the atom has been seen before - COULD LEAD TO EXTRA CHAINS BEING ADDED THAT ARE NOT USED.
                            $fpr[$addchain]=$bond_list_ptr[1][$i][$j2];     # ptr of first atom in new chain
                            $cord[$addchain]=$addchain;                     # chain value/number
                                                                            # add chain tells you how many chains that need to be prioritised.
                           # printf "Atom i: %i branch: %i from: %i chain id: %i\n", $lmat[1][$fpr[$addchain]], $addchain, $lmat[1][$i], $lmat[8][$fpr[$addchain]];
                            #printf "Atom i: %i branch: %i from: %i chain id: %i\n", $dist_ca[$fpr[$addchain]], $dist_ca[$i], $lmat[1][$i];
                            }
                        }
                    }
                }
#    printf "b %i %i %i\n", $lmat[1][$i], $lmat[8][$i], $addchain;

    #       i have list of atoms that need to be prioritised
        if($addchain > 1)
        {
           # printf "BRANCH branch: %i from: %i\n", $addchain, $lmat[1][$i], ;
            $nors=0;
            sub_priority();

############## ABOVE MATRIX CREATED
#############       END ASSIGNING PRIORITY  ##########
    
        }
    
#    printf "%i %i %i\n", $lmat[1][$i], $lmat[8][$i], $i;

#    $lmat[8][$i]=$t_chains;      #ca chain will be chain 1
#    for ($j=$i+1 ; $j < $nhat ; $j++)
    
###     The bit below takes the information about priorities coming out of the priority subroutine and uses this to increment chains.
#       First bit just establishes which branching atoms are to be considered [this information is in addchains, but would not be present for non-branch points]
#######################################################################################################
#   If the atom
#    for all atoms after CA -
#        check all connnected atoms 1 bond away - 1 at a time
#            if an atom 1 bond away from a given atom is this atom (i.e. it is the next atom along the chain)
#                if that atom has not been seen - or it has been seen and appears after this atom
#                    then do something unless
#                        the atom already belongs to a chain that has higher priority than this chain
#    tm_n = 1 then does the something.
#
#######################################################################################################
    for ($j=2 ; $j < $nhat ; $j++)
    {
        #printf "%i %i %i\n", $lmat[1][$j], $lmat[8][$j], $dist_ca[$j];
        $m_tn = 0;
        for ($j2=1 ; $j2 <= $bond_lst_l[1][$j] ; $j2++)         # check how many branching atoms > $i - > if less than $i then not branch
        {
            if ($bond_list[1][$j][$j2] == $lmat[1][$i] )                #
            {
                    # if atoms are connected and the atom has not been seen before - or if it has been seen if it appears after this atom then print
                if (($lmat[11][$j] == 0 ) || ($lmat[11][$j] > 0 && $lmat[11][$j] > $lmat[11][$i]) )
                {
                    $m_tn=1;
                #printf "connected2 tot atoms: %i atom: %i one-bond to: %i - order of j, %i order of i %i addchain %i\n", $bond_lst_l[1][$j] ,$bond_list[1][$j][$j2], $lmat[1][$j], $lmat[8][$j], $lmat[8][$i], $addchain;
                    if ($lmat[11][$j] > 0)    # atom has been seen before
                    {
                        if ($lmat[8][$j] < $lmat[8][$i])    # if the priority of that previously seen atom is higher than current chain then ignore
                        {
                            $m_tn=0;
                 #         printf "bottom of ring? %i %i %i %i\n", $lmat[1][$j] , $lmat[1][$i], $lmat[8][$j] , $lmat[8][$i];
                        }
                    }
                }
            }
        }
#        printf "x %i %i %i\n", $lmat[1][$i], $lmat[8][$i], $m_tn;
#######################################################################################################
# add that atom to this chain [i.e. give it the same chain ID]
# and give it a distance to ca - which tells the above that it has been seen now
#    in cases where you have more than 1 chain being evaluated
#        check each chain and only apply this logic to the chain with the highest priority
#            the other chains with lower priority will get an incremented chain ID.
#
#######################################################################################################

        if ($m_tn == 1 )
        {
#   don't update the chain id if the chain id of that atom is lower than the chain id of this atom.
            $assigned=0;
            if ($addchain <= 1)
            {
                $lmat[8][$j]=$lmat[8][$i];      #same chain id as what it is connected to.
                $assigned=1;
            }
                $lmat[11][$j]=$lmat[11][$i]+1;  # this should be moot
                #printf "connected j %i, i %i i %i tn %i ord j %i i %i\n", $lmat[1][$j], $lmat[1][$i], $i, $m_tn, $lmat[8][$j], $lmat[8][$i];
            if ($addchain >= 1)
            {

                $c_chains=$addchain-1;
                #$c_skip=0;
                #printf "j: %i \n", $j;

                for ($p = 1 ; $p <= $addchain ; $p++)
                {
                    #printf "addchain  %i %i, j %i, i %i, cord[, %i, fpr, %i\n", $p, $lmat[1][$fpr[$p]], $lmat[8][$j], $lmat[8][$i], $cord[$p], $fpr[$p] ;

                    if ($fpr[$p] == $j)
                    {
                        if ($cord[$p] == 1)                                 #   There is an issue here that means that if the second priority atom is the continuation of the chain it will not be done correctly, and instead it is given a new chain ID. This is OK as there is a bit later that fixes this. But the logic is not perfect here. You would have to keep track of if $lmat[8][$i] has been given to anyone and then only make a new chain if it has.
                        {
                            $lmat[8][$j]=$lmat[8][$i];
                            $assigned=1;

                        }else{
                     #       if (($lmat[8][$j] > 0) )
                     #       {
                     #         printf "skip with %i\n", $lmat[8][$j];
                     #       }
                            if ($lmat[8][$j] == 0)
                            {
                     #           printf "? %i %i %i\n", $lmat[8][$j], $t_chains, $cord[$p]-1;
                                $lmat[8][$j] = $t_chains+$cord[$p]-1;           #   THIS WILL GENERATE SOME EMPTY CHAINS WHICH WILL GET FIXED LATER
                            }
                        }
                            
                    #printf "Branch: bpoint %i branch %i to atom: %i with chain id %i total chains %i cord p %i p %i %i\n", $lmat[1][$i], $lmat[8][$i], $lmat[1][$j], $lmat[8][$j], $t_chains, $cord[$p], $p, $assigned;
                   # printf "p: %i, swapped with cordp %i \n", $p, $cord[$p];

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
#for ($i=1 ; $i < $nhat ; $i++)
#{
#    printf "dists: %i %i %i %i\n", $dist_ca[$i], $lmat[1][$i], $lmat[11][$i], $lmat[8][$i];
#}
#
#printf "$at_ca\n";
#for ($j=1 ; $j <= $t_chains ; $j++)         ######### f
#{
#    printf "CHAIN %i (%i): ", $j, $cxl[$j];
#    for ($k = 0 ; $k < $cxl[$j] ; $k++)
#    {
#        printf "%i ", $cx[$j][$k], $cxl[$j];
#    }
#    printf "\n";
#
#}
############################################
############################################  JOIN CHAINS THAT ARE BROKEN

$nj=0;
for ($j=1 ; $j <= $t_chains ; $j++)
{
    #printf "CHAIN %i: ", $j;
    for ($k = 0 ; $k < $cxl[$j] ; $k++)
    {
        if ($k == $cxl[$j]-1)
        {
            for ($j2=1 ; $j2 <= $t_chains ; $j2++)
            {
                for ($k2 = 0 ; $k2 < $cxl[$j2] ; $k2++)
                {
                    if ($k2 == 0 && $j != $j2)
                    {
                        ####################        LAST: $cx[$j][$k], FIRST(s): $cx[$j2][$k2]      ############

                        #printf "Compare: LAST %i %i TO FIRST %i %i\n", $cx[$j][$k], $lmat[11][$cx_1[$j][$k]], $cx[$j2][$k2], $lmat[11][$cx_1[$j2][$k2]] ; #, $cx_1[$j][$k];
                        if ($lmat[11][$cx_1[$j][$k]]+1 == $lmat[11][$cx_1[$j2][$k2]])
                        {
                            for ($i=1 ; $i < $nhat ; $i++)
                            {
                                if ($lmat[1][$i] == $cx[$j][$k])        # LAST ATOM CHECK ALL 1 BOND CONNECTIONS
                                {
                                    for ($j3=1 ; $j3 <= $bond_lst_l[1][$i] ; $j3++)       # how many entries in x-bond list - go through each
                                    {
                                        if ($bond_list[1][$i][$j3] == $cx[$j2][$k2] && $bond_list[1][$i][$j3] < 1000)
                                        {
                                            #printf "\n1 1 bond OF last %i  to first of other chain %i\n", $bond_list[1][$i][$j3], $cx[$j2][$k2];
                                            #printf "\n join %i to end of %i \n", $j2, $j;
                                            $bj[$nj]=$j;
                                            $ej[$nj]=$j2;
                                            $nj++;
                                        }
                                        
                                        
                                    }
                                    #
                                }
                            }
                        }

                        ####################
                        #printf "I'm first %i ", $cx[$j2][$k2]; #, $cx_1[$j][$k];
                    }

                }
            }
        }
#        printf "%i ", $cx[$j][$k]; #, $cx_1[$j][$k];
    }
#    printf "\n";
}
#
for ($j=1 ; $j <= $t_chains ; $j++)
{
    $gc[$j]=0;
}

for ($my_j=1 ; $my_j <= $t_chains ; $my_j++)         #########
{
    $addchain=0;
    for ($i = 0; $i <$nj ; $i++)
    {
        #printf "A %i\n", $addchain;
        if ($bj[$i] == $my_j)
        {
            $lastpt = $cxl[$my_j]-1;
            
         #   printf "JOIN THIS CHAIN %i (at: %i -adchain %i): \n", $my_j, $cx[$my_j][$lastpt], $addchain;
            $addchain++;
            $ptr1 = $cx_1[$ej[$i]][0];
            $fpr[$addchain]=$ptr1;
            $cord[$addchain]=$ej[$i];       # chain id of the attached thing
        }
        
    }
    if ($addchain > 1)
    {
        $nors=1;
            sub_priority();
        #printf "BRANCH branch: %i from: %i\n", $addchain, $lmat[1][$i], ;

        for ($i = 0; $i <$nj ; $i++)
        {
            if ($bj[$i] == $my_j)
            {
                
                #printf "JOIN THIS CHAIN %i (at: %i -adchain %i): \n", $my_j,$ej[$i], $addchain;
                #            $cxl[$j]=0;
                #            $cxp[$j]=$cxp[$ej[$i]];
                #            $cxptr[$j]=$cxptr[$ej[$i]];
                for ($ai = 2 ; $ai <= $addchain ; $ai++)
                {
                    if ($cord[$ai] == $ej[$i])
                    {
                        $bj[$i]=-1;
                    }
                 #   printf "Y %i %i\n", $cord[1], $addchain;
                }
            }
        }
        
    }
    #printf "\n";
    
}




for ($j=1 ; $j <= $t_chains ; $j++)         #########
{

    for ($i = 0; $i <$nj ; $i++)
    {
        if ($bj[$i] == $j)
        {
#            printf "JOIN THIS CHAIN %i: \n", $j;
#            $cxl[$j]=0;
            $cxp[$j]=$cxp[$ej[$i]];
            $cxptr[$j]=$cxptr[$ej[$i]];
            for ($k2 = 0 ; $k2 < $cxl[$ej[$i]] ; $k2++)
            {
 #               printf "%i ", $cx[$ej[$i]][$k2]; #, $cx_1[$j][$k];
                $l=$cxl[$j];
                $cx[$j][$l]=$cx[$ej[$i]][$k2];
                $cx_1[$j][$l]=$cx_1[$ej[$i]][$k2];
                $cxl[$j]++;
            }
            $gc[$ej[$i]]=-1;
        }

    }
#    printf "\n";
    
}


#
for ($j=1 ; $j <= $t_chains-1 ; $j++)         ######### f
{
#    printf "xj %i \n", $j;
    for ($j2=$j+1 ; $j2 <= $t_chains ; $j2++)         ######### for
    {
#        printf "xj2 %i \n", $j2;

        if ($gc[$j] < $gc[$j2] )  # overwrite this entry with last and truncate
        {
            for ($k = 0 ; $k < $cxl[$j2] ; $k++)
            {
                $cx[$j][$k]=$cx[$j2][$k]; #, $cx_1[$j][$k];
                $cx_1[$j][$k]=$cx_1[$j2][$k];
            }
            $cxl[$j] = $cxl[$j2] ;
            $cxp[$j]=$cxp[$j2];
            $cxptr[$j]=$cxptr[$j2];
            $tgc = $gc[$j];
            $gc[$j]=$gc[$j2]    ;
            $gc[$j2] = $tgc;
        }
    }
}
#for ($j=1 ; $j <= $t_chains ; $j++)
#{
#    printf "xsCHAINss %i: (%i)", $j, $gc[$j];
#    for ($k = 0 ; $k < $cxl[$j] ; $k++)
#    {
#        printf "%i ", $cx[$j][$k]; #, $cx_1[$j][$k];
#    }
#    printf "\n";
#
#}


$e_chains = 0;
for ($j=1 ; $j <= $t_chains ; $j++)         ######### for
{
    if ($gc[$j] == -1)
    {
        $e_chains = $j;
        last;
    }
}
if ($e_chains > 0) {$t_chains = $e_chains-1};
#for ($j=1 ; $j <= $t_chains ; $j++)         ######### f
#{
#    printf "CHAIN %i (%i): ", $j, $cxl[$j];
#    for ($k = 0 ; $k < $cxl[$j] ; $k++)
#    {
#        printf "%i ", $cx[$j][$k], $cx_1[$j][$k];
#        #printf "%i ", $cx[$j][$k];
#    }
#    printf "\n";
#
#}



############################################
############################################ END JOINGING CHAINS
######### swap chains based on order
for ($cxj=1 ; $cxj < $t_chains ; $cxj++)
{
    $addchain = 1;
    $ptr1 = $cx_1[$cxj][0];         # first atom in chain XJ

    $fpr[$addchain]=$ptr1;
    $cord[$addchain]=$cxj;
    for ($cxj2=$cxj+1 ; $cxj2 <= $t_chains ; $cxj2++)
    {
        if ($cxl[$cxj] > 0 && $cxl[$cxj2] > 0)
        {
        $ptr2 = $cx_1[$cxj2][0];    # first atom in chain XJ2
        if ($lmat[11][$ptr1] == $lmat[11][$ptr2])
        {
            $addchain++;
            $fpr[$addchain]=$ptr2;
            $cord[$addchain]=$cxj2;
            #printf "Compare: FIRST %i TO FIRST %i. %i %i %i %i\n", $lmat[1][$ptr1], $lmat[1][$ptr2], $lmat[11][$ptr2], $lmat[11][$ptr1], $cxj, $cxj2 ; #, $cx_1[$j][$k];
        }
        }
    }
############
            #############
            #############
            if($addchain > 1)
            {
#                for ($p = 1 ; $p <= $addchain ; $p++)       # PRINTING ONLY
#                {
#                    for ($cxj2=1 ; $cxj2 <= $t_chains ; $cxj2++)
#                    {
#                        if ($fpr[$p] == $cx_1[$cxj2][0])    # first atom in that chain matches fpr
#                        {
#                            printf "chain %i has %i priority\n", $cxj2, $p;
#
#                        }
#                    }
#
#                }
                $nors=1;
                $ed=0;
                sub_priority();
                ######################      if two chains are indentical then order these based on what they are connected to   #####################
                if ($ed==0)     # they are otherwise identical
                {
                    for ($p = 1 ; $p < $addchain ; $p++)       # PRINTING ONLY
                      {
                          for ($cxj2=1 ; $cxj2 <= $t_chains ; $cxj2++)
                          {
                              if ($fpr[$p] == $cx_1[$cxj2][0])    # first atom in that chain matches fpr
                              {
                                 $tx1=$lmat[9][$fpr[$p]];
                                  for ($ix=1 ; $ix < $nhat ; $ix++)
                                  {
                                      if ($lmat[1][$ix] == $tx1)
                                      {
                                         # printf "prioirty of attached atom is: %i\n", $lmat[8][$ix];
                                          $p1pr = $lmat[8][$ix];
                                          $px1=$p;
                                      }
                                  }
                                 # printf "chain %i has %i priority\n", $cxj2, $p;

                              }
                          }
                          for ($p2 = $p+1 ; $p2 <= $addchain ; $p2++)       # PRINTING ONLY
                            {
                                for ($cxj2=1 ; $cxj2 <= $t_chains ; $cxj2++)
                                {
                                    if ($fpr[$p2] == $cx_1[$cxj2][0])    # first atom in that chain matches fpr
                                    {
                                       $tx2=$lmat[9][$fpr[$p2]];
                                        for ($ix=1 ; $ix < $nhat ; $ix++)
                                        {
                                            if ($lmat[1][$ix] == $tx2)
                                            {
                                        #        printf "prioirty of attached atom is: %i\n", $lmat[8][$ix];
                                                $p2pr = $lmat[8][$ix];
                                                $px2=$p2;
                                            }
                                        }
                                       # printf "chain %i has %i priority\n", $cxj2, $p2;

                                    }
                                }
                            }
                          if ($p1pr > $p2pr)
                          {
                              #printf "compared %i with %i and chain priority is: %i %i %i\n", $px1, $px2, $p1pr, $p2pr, $addchain;
                              ($fpr[$px1], $fpr[$px2]) = ($fpr[$px2], $fpr[$px1]);
                          }
                      }
                }
####   END                ######################      if two chains are indentical then order these based on what they are connected to   #####################

                #printf "ed %i\n", $ed;
                for ($p = 1 ; $p <= $addchain ; $p++)
                {
                    for ($cxj2=1 ; $cxj2 <= $t_chains ; $cxj2++)
                    {
                        if ($fpr[$p] == $cx_1[$cxj2][0])    # first atom in that chain matches fpr
                        {
                            $i=$cxj2;
                            #printf "chain %i has %i priority. length %i. cord: %i\n", $cxj2, $p, $cxl[$i], $cord[$p];
                            if ($cxj2 != $cord[$p])
                            {
                             #   printf "Swap %i %i\n", $cord[$p], $cxj2;
                                $i=$cxj2;
                                $j=$cord[$p];
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
                                         $cx_1[$i][$k]     =  $cx_1[$j][$k];
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
                }
            #printf "BRANCH branch: %i from: %i\n", $addchain, $lmat[1][$i], ;


        
            }
    
}
#for ($j=1 ; $j <= $t_chains ; $j++)         ######### f
#{
#    printf "A1: CHAIN %i: ", $j;
#    for ($k = 0 ; $k < $cxl[$j] ; $k++)
#    {
#        printf "%i ", $cx[$j][$k];
#    }
#    printf "\n";
#
#}

##   remove empty chains
for ($i = 1 ; $i < $t_chains ; $i++)
{
    for ($j = $i+1 ; $j <= $t_chains ; $j++)
    {
        if ($cxl[$i] == 0)
#            if ($cxl[$j] > $cxl[$i])
        {
#            printf "Warning orphan chain present.\n";
            $tvl1       = $cxl[$i];
            $tprt       = $cxp[$i];
            for ($k = 0 ; $k < $cxl[$i] ; $k++)
            {
                $tvcx1[$k]  =   $cx[$i][$k];
                $tvcx5[$k]  = $cx_1[$i][$k];
            }
            for ($k = 0 ; $k < $cxl[$j] ; $k++)
            {
                   $cx[$i][$k]     =    $cx[$j][$k];
                 $cx_1[$i][$k]     =  $cx_1[$j][$k];
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
$shrink=0;

for ($j=1 ; $j <= $t_chains ; $j++)         ######### f
{
    if ($cxl[$j] == 0)
    {
            $shrink++;
    }
}
$t_chains = $t_chains-$shrink;
#for ($j=1 ; $j <= $t_chains ; $j++)         ######### f
#{
#    printf "A1: CHAIN %i: ", $j;
#    for ($k = 0 ; $k < $cxl[$j] ; $k++)
#    {
#        printf "%i ", $cx[$j][$k];
#    }
#    printf "\n";
#
#}
###############   UPDATE LMAT[8]  ###########
for ($j=1 ; $j <= $t_chains ; $j++)
{
      for ($k = 0 ; $k < $cxl[$j] ; $k++)
       {
           for ($i=1 ; $i < $nhat ; $i++)
           {
               if ($lmat[1][$i] == $cx[$j][$k])
               {
                    $lmat[8][$i] = $j ;
               }
           }
       }
}
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
#           printf "%i %i %i ", $cx[$j][$k], $cx_1[$j][$k];
            $p_cx[$j][$t2]=$cx[$j][$k];
    }
#    printf "\n";
}
################ COMPARE POSITIONALLY
for ($x=0 ; $x <= $maxlmat ; $x++)
{
    $count_matl=0;
    for ($j=1 ; $j <= $t_chains ; $j++)
    {
            #printf "%i\t", $p_cx[$j][$x];
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
           # printf "%s \n %s\n", $fil_c[$i], $fields[1];

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
         $allnh++;
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
            #printf "3at no: %i carb: %i %i \n", $lmat[1][$k], $hmat[2][$i], $hmat[4][$i];
            $hmat[20][$i]=$lmat[8][$k];          #same chain, used later for printing
            if ($lmat[8][$k] == 1 && $c_ch_at[$lmat[11][$k]] == 1)
            {
                    if ($hmat[4][$i] == 1)
                    {
                        $hmat[6][$i] = sprintf("H%s  ", $gnam[$lmat[11][$k]]);
                        #printf "3at no: %i carb: %i %i \n", $lmat[1][$k], $hmat[2][$i], $hmat[4][$i];
                        #printf "%s %i\n", $hmat[6][$i], $hmat[1][$i];
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
                                    
                            for ($xi=$inat ; $xi <=$endat ; $xi++)         #   print hydrogens
                            {
                                @fields = split (' ', $fil_c[$xi]);
                                if ($xi > $atstart)
                                {
                             #       printf "%f %f %f\n",  $fields[2], $fields[3], $fields[4];
                             #       printf "c: %i %i %i %i %i %i\n",  $fields[0], $lmat[1][$k], $lmat[9][$i], $at_n, $lmat[30][$k],  $hmat[1][$i];

                                    if($fields[0] == $lmat[1][$k])        # PUT SUBSTIUENT WHERE THE CENTRAL ATOM IS AND GIVE LOWEST PRIORITY
                                    {
                                        @at3 = ($fields[2], $fields[3], $fields[4]);
                                     #   printf "inside %f %f %f\n",  $fields[2], $fields[3], $fields[4];

                                    }
                                    if($fields[0] == $at_n && $lmat[9][$k] == 0)    # first nitrogen you hit - shoudl be 1.
                                    {
                                        @at2 = ($fields[2], $fields[3], $fields[4]);
                                    }
                                    
                                    if($fields[0] == $lmat[9][$k] && $lmat[9][$k] > 0)
                                    {
                                        #print $lmat[9][$i] $fields0;
                                        @at2 = ($fields[2], $fields[3], $fields[4]);
                                    }
                                    if($fields[0] == $lmat[30][$k])        # first hydrogen atom encountered on that carbon
                                    {
                                        @at1 = ($fields[2], $fields[3], $fields[4]);
                                    }
                                    if($fields[0] == $hmat[1][$i])
                                    {
                                        @at4 = ($fields[2], $fields[3], $fields[4]);
                                    }
                                }
                            }
                            #print "@at1, @at2, @at3, @at4\n";
                            &getrs(@at1, @at2, @at3, @at4, $det);
                            if ($det < 1)
                            {
                             #   printf "Clockwise.\n";
                            }
                            if ($det > 1)       # should only be here on second atom of 2
                            {
                              #  printf "Anti-clockwise. \n";
                                $anti=1;
                            }
                            
                            
                        }
                        ##################    END   for CIP   ##################

                        $hmat[6][$i] = sprintf("H%s%i ", $gnam[$lmat[11][$k]], $hmat[5][$i]+1);
                        ##################     for CIP2   ###################
                        if ($hmat[5][$i] < $hmat[4][$i])          # if first of a pair i store atom names for later swap
                        {
                            $lmat[31][$k]=$hmat[6][$i];
                            $lmat[32][$k]=$i;
                        }
                        if ($anti == 1 && $hmat[5][$i] == $hmat[4][$i])          # last of a pair, if this atom being considered is clockwise (means HB3 is in position of HB2), then overwrite
                        {                                                        # if the last of a pair is the pro-R then swap
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
                         #printf "central: %i back: %i  h1: %i h2:%i\n", $lmat[1][$k], $lmat[9][$k], $hmat[1][$i], $lmat[30][$k];
                         
                         for ($xi=$inat ; $xi <=$endat ; $xi++)         #   print hydrogens
                         {
                             @fields = split (' ', $fil_c[$xi]);
                                if ($xi > $atstart)
                                 {
                                     #printf "%f %f %f\n",  $fields[2], $fields[3], $fields[4];
                                     #printf "c: %i %i %i %i %i %i\n",  $fields[0], $lmat[1][$k], $lmat[9][$i], $at_n, $lmat[30][$k],  $hmat[1][$i];

                                     if($fields[0] == $lmat[1][$k])        # PUT SUBSTIUENT WHERE THE CENTRAL ATOM IS AND GIVE LOWEST PRIORITY
                                     {
                                         @at3 = ($fields[2], $fields[3], $fields[4]);
                                         #printf "inside %f %f %f\n",  $fields[2], $fields[3], $fields[4];

                                     }
                                     if($fields[0] == $at_n && $lmat[9][$k] == 0)    # first nitrogen you hit - shoudl be 1.
                                     {
                                         @at2 = ($fields[2], $fields[3], $fields[4]);
                                     }
                                     
                                     if($fields[0] == $lmat[9][$k] && $lmat[9][$k] > 0)
                                     {
                                         #print $lmat[9][$i] $fields0;
                                         @at2 = ($fields[2], $fields[3], $fields[4]);
                                     }
                                     if($fields[0] == $lmat[30][$k])        # first hydrogen atom encountered on that carbon
                                     {
                                         @at1 = ($fields[2], $fields[3], $fields[4]);
                                     }
                                     if($fields[0] == $hmat[1][$i])
                                     {
                                         @at4 = ($fields[2], $fields[3], $fields[4]);
                                     }
                                 }
                             }
                             #print "@at1, \n@at2, \n@at3, \n@at4\n";
                             &getrs(@at1, @at2, @at3, @at4, $det);
                             if ($det < 1)
                             {
                              #   printf "Clockwise. 2\n";
                             }
                             if ($det > 1)       # should only be here on second atom of 2
                             {
                               #  printf "Anti-clockwise. 2\n";
                                 $anti=1;
                             }
                         
                         
                     }
                     ##################    END   for CIP   ##################
                    $hmat[6][$i] = sprintf("H%s%i%i", $gnam[$lmat[11][$k]], $lmat[8][$k], $hmat[5][$i]+1);
                    ##################     for CIP2   ###################
                    if ($hmat[5][$i] < $hmat[4][$i])          # if first of a pair i store atom names for later swap
                    {
                        $lmat[31][$k]=$hmat[6][$i];
                        $lmat[32][$k]=$i;
                    }
                    if ($anti == 1 && $hmat[5][$i] == $hmat[4][$i])          # last of a pair, if this atom being considered is clockwise (means HB3 is in position of HB2), then overwrite
                    {                                                        # if the last of a pair is the pro-R then swap
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

#    if ($hmat[2][$i] == 3 && $nterm == 0)
#    {
#        $hmat[6][$i] = sprintf("H   ");     #special case for amide H
#        $hmat[20][$i]=1;
#    }
    if ($hmat[2][$i] == 1 && $nterm == 1)
    {
        if ($allnh == 1)
        {
        $hmat[6][$i] = sprintf("H   ");     #special case for amide H
        }else{
        $hmat[6][$i] = sprintf("H%i  ", $hmat[21][$i]);     #special case for amide H
        }
        $hmat[20][$i]=1;
    }
        #printf "n %s %i %i\n", $hmat[6][$i], $hmat[1][$i], $hmat[2][$i];

}
for ($i=1 ; $i <=$hhat ; $i++)        # get pseudoatoms
{
    if ($hmat[5][$i] > 0)
    {
            printf "%i %s\n", $hmat[1][$i], $hmat[6][$i];
    }
}
##############   20220614 - for re-ordering.
$k=0;
for ($i=1 ; $i < $nhat ; $i++)
{
    $k++;
    $orgmat[$k]=$lmat[1][$i];
    #printf "R: %i %s\n", $lmat[1][$i], $lmat[13][$i];
    for ($j=1 ; $j <=$hhat ; $j++)        # get pseudoatoms
    {
        if ($hmat[5][$j] > 0 && $hmat[2][$j] == $lmat[1][$i])
        {
                $k++;
                $orgmat[$k]=$hmat[1][$j];
              #  printf "R: %i %s\n", $hmat[1][$j], $hmat[6][$j];
        }
    }

}
$maxk=$k;
my @sorted_orgmat = sort { $a <=> $b } @orgmat;
for ($k=1 ; $k<=$maxk ; $k++)
{
 #   printf "%i %i\n", $orgmat[$k], $sorted_orgmat[$k];
}
$k=0;
for ($i=1 ; $i < $nhat ; $i++)
{
    $k++;
    printf "S: %i %s %i\n", $lmat[1][$i], $lmat[13][$i], $sorted_orgmat[$k];
    for ($j=1 ; $j <=$hhat ; $j++)        # get pseudoatoms
    {
        if ($hmat[5][$j] > 0 && $hmat[2][$j] == $lmat[1][$i])
        {
                $k++;
                printf "S: %i %s %i\n", $hmat[1][$j], $hmat[6][$j], $sorted_orgmat[$k];
        }
    }

}
############### PRINT #############
$i=0;
for ($xi=$inat ; $xi <=$endat ; $xi++)        #   print hydrogens
{
    @fields = split (' ', $fil_c[$xi]);
    if ($xi > $atstart)
    {
        # LOOKING FOR STANDARD N/C TERMINAL ATOMS
        if ($fields[0] == $at_tco)
        {
            #print "found co\n";
            $bbmat[1][$i] = $fields[0];
            $bbmat[2][$i] = sprintf("O2");
            $i++;
        }
        if ($fields[0] == $at_hn2)
        {
           # print "found HN2\n";
            $bbmat[1][$i] = $fields[0];
            $bbmat[2][$i] = sprintf("H2");
            $i++;

        }
        if ($fields[0] == $at_hn3)
        {
            print "found H3\n";
            $bbmat[1][$i] = $fields[0];
            $bbmat[2][$i] = sprintf("H3");
            $i++;
        }
        # LOOKING FOR BACKBONE HN N C' and (C')OTERMINAL ATOMS
        if ($fields[0] == $at_hn && $at_hn2 > 0)
        {
            #print "found H\n";
            $bbmat[1][$i] = $fields[0];
            $bbmat[2][$i] = sprintf("H1");
            $i++;
        }elsif($fields[0] == $at_hn && $at_hn2 < 0)
        {
            $bbmat[1][$i] = $fields[0];
            $bbmat[2][$i] = sprintf("HN");
            $i++;
        }
            
        if ($fields[0] == $at_c)
        {
            #print "found C'\n";
            $bbmat[1][$i] = $fields[0];
            $bbmat[2][$i] = sprintf("C");
            $i++;
        }
        if ($fields[0] == $at_n)
        {
            #print "found N\n";
            $bbmat[1][$i] = $fields[0];
            $bbmat[2][$i] = sprintf("N");
            $i++;
        }
        if ($fields[0] == $at_co)
        {
            #print "found CO\n";
            $bbmat[1][$i] = $fields[0];
            $bbmat[2][$i] = sprintf("O");
            $i++;
        }


    }
}
$allbb=$i;
#print ("$allbb\n");

#$mol2 = sprintf("block_%s", $inputfile);

open(OUT, '>', $outputfile) || die "Can't open output file";
#        printf OUT "%s",$head[$i];

$endb=0;
for ($j = 1 ; $j <= $numlin ;  $j++)
{
    @fields = split (' ', $fil_c[$j]);

     if ($fil2_c[$j] == 3)
     {
            $endb=0;
            for ($i=0 ; $i < $allbb ; $i++)      ### just for printing
            {
             if ($bbmat[1][$i] == $fields[0])
             {
                # print ("$bbmat[2][$i]\n");
                 $bbmat[2][$i]=sprintf("%4s", $bbmat[2][$i]);
                 substr($fil_c[$j],8,4)=$bbmat[2][$i];
                 printf OUT "%s\n", $fil_c[$j];
                 $endb=1
             }
            }

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
            if ($j == $numlin )
            {
                printf OUT "%s", $fil_c[$j];
            }else{
             printf OUT "%s\n", $fil_c[$j];
            }
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


sub sub_priority
{
    #printf "insub %i %i\n", $lmat[1][1], $at_ca;
    
    ################ Check priority of hits and sort based on hit #####################
    ###########      Check priority of atom itself.                 ######## needs to be done separately as it is not part of the bond_list.
    
    for ($p = 1 ; $p <= $addchain ; $p++)
    {
        $x=0;
        $i2 = $fpr[$p];
        $pr_mat[$p][$x] =  $lmat[10][$i2];
        for ($x = 1; $x < $nhat ; $x++)         #### SORTS THE PRIORITIES AT EACH BOND DISTANCE SO THAT THE ATOM WITH THE HIGHEST PRIORITY IS ON TOP OF PR_MAT
        {
            $pr_mat[$p][$x] =  0;
            if ($bond_lst_l[$x][$i2] > 0)
            {
                
                $i3             = $bond_list_ptr[$x][$i2][1];
                $maxp           = $lmat[10][$i3];                              # the first entry
                $pr_mat[$p][$x] = $lmat[10][$i3];
                
                for ($j=1 ; $j <= $bond_lst_l[$x][$i2] ; $j++)       # how many entries in x-bond list - go through each
                {
                    $i3=$bond_list_ptr[$x][$i2][$j];
                    if ($lmat[10][$i3] > $maxp )
                    {
                        $pr_mat[$p][$x] = $lmat[10][$i3];
                        $maxp           = $lmat[10][$i3];
                    }
                }
                # printf "step 1: Atom i: %i Bonds away: %i Priority: %i\n", $lmat[1][$fpr[$p]], $x, $maxp, ;
                
            }
        }
        $fpr_i[$i]=0;
    }
    #                                            for ($p = 1 ; $p <= $addchain ; $p++)
    #                                            {
    #                        #                            printf "Atom i: %i Bonds away: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$p]], $x, $pr_mat[$p][$x], $p;
    #                                                    printf "stuff3 A i: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$p]], $cord[$p], $p;
    #                                                    #if ($pr_mat[$p][$x] == 0) {printf "should exit here saying that assignment is arbitrary\n"};
    #                                            }
    for ($p = 1 ; $p <= $addchain-1 ; $p++)
    {
        for ($p2 = 1+$p ; $p2 <= $addchain ; $p2++)
        {
            ###     STEP 1: ##########  FIRST DETERMINE PRIORITY OF SELF ###############    IF PRIORITY OF SELF IS HIGHER/LOWER THEN ORDER CHAINS AND MOVE TO NEXT CHAIN
            #printf "Self Atom i: %i %i Priority: %i chain: %i\n", $lmat[1][$fpr[$p]],$lmat[1][$fpr[$p2]], $pr_mat[$p][0], $pr_mat[$p2][0];
            
            if ($pr_mat[$p2][0] != $pr_mat[$p][0])
            {
                #    printf "swap1\n";
                if ($pr_mat[$p2][0] > $pr_mat[$p][0])
                {
                    #    printf "swap2\n";
                    # $temp_ch            = $cord[$p];
                    # $cord[$p]           = $cord[$p2];
                    # $cord[$p2]          = $temp_ch;
                    $temp_fp            = $fpr[$p];
                    $fpr[$p]            = $fpr[$p2];
                    $fpr[$p2]           = $temp_fp;
                    
                    for ($x2 = 0; $x2 < $nhat ; $x2++)
                    {
                        $temp_val           = $pr_mat[$p][$x2];
                        $pr_mat[$p][$x2]    = $pr_mat[$p2][$x2];
                        $pr_mat[$p2][$x2]   = $temp_val;
                    }
                }
                $ed=1;
                next;
            }
            ###     STEP 2: ########## IF SELF CANNOT DISTINGUISH PRIORITY  EVALUATE ALL CONNECTED ATOMS AT EACH POSITION ##############
            for ($x = 1; $x < $nhat ; $x++)
            {
                #printf "Atom i: %i Bonds away: %i Priority: %i chain: %i, longest:L %i\n", $lmat[1][$fpr[$p]], $x, $pr_mat[$p][$x], $p, $longest;
                #printf "CF A i: %i Bonds away: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$p2]], $x, $pr_mat[$p2][$x], $p2;
                #if ($pr_mat[$p][$x] == 0) {printf "should exit here saying that assignment is arbitrary\n"};
                $i2=$fpr[$p];
                $i22=$fpr[$p2];
                $fine=0;
                $longest = $bond_lst_l[$x][$i2];
                if ($bond_lst_l[$x][$i2] < $bond_lst_l[$x][$i22] ) {$longest = $bond_lst_l[$x][$i22]};
                #printf "ll %i\n",$longest;
                
                ###     SECOND DETERMINE PRIORTY TO ALL CONNECTED ATOMS
                
                if ( ($pr_mat[$p2][$x] == $pr_mat[$p][$x]) && ($longest > 1))        ####### IF THE ATOM WITH THE HIGHEST PRIORITY X BONDS AWAY IS NOT DIFFERENT
                {
                    #printf "Current mess.\n";
                    # MAKE 2 DUMMY LISTS - DL1, DL2, INITIALISE LISTS TO SAME LENGTH
                    $dll = $bond_lst_l[$x][$i22];
                    if ( $bond_lst_l[$x][$i2] > $bond_lst_l[$x][$i22] ){ $dll = $bond_lst_l[$x][$i2]};
                    
                    for ($dj = 0; $dj <= $nhat ; $dj++)     # sort will reorder everythign in this array, so it must be 0 for all possible length - alternatively if the sort function is truncated it would also work
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
                    for ($dj = 0; $dj < $dll ; $dj++)
                    {
                        #printf "Comparing %i %i.\n", $dl1s[$dj] , $dl2s[$dj];
                        
                        if ($dl1s[$dj] > $dl2s[$dj])
                        {
                            #   printf "i2 %i %i wins over %i. \n", $lmat[1][$i2], $dl1s[$dj] , $dl2s[$dj];
                            $ed=1;
                            $fine=1;
                            last;
                            
                        }
                        if ($dl1s[$dj] < $dl2s[$dj])
                        {
                            
                            #  printf "i22 %i %i wins over %i. %i %i\n", $lmat[1][$i22], $dl2s[$dj] , $dl1s[$dj], $p2, $p;
                            $ed=1;
                            
                            # $temp_ch            = $cord[$p];
                            # $cord[$p]           = $cord[$p2];
                            # $cord[$p2]          = $temp_ch;
                            $temp_fp            = $fpr[$p];
                            $fpr[$p]            = $fpr[$p2];
                            $fpr[$p2]           = $temp_fp;
                            
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
                    
                }
                
                if($pr_mat[$p2][$x] != $pr_mat[$p][$x])
                {
                    $ed=1;
                    
                    #printf "Second. \n";
                    if ($pr_mat[$p2][$x] > $pr_mat[$p][$x])
                    {
                        # $temp_ch            = $cord[$p];
                        # $cord[$p]           = $cord[$p2];
                        # $cord[$p2]          = $temp_ch;
                        $temp_fp            = $fpr[$p];
                        $fpr[$p]            = $fpr[$p2];
                        $fpr[$p2]           = $temp_fp;
                        
                        for ($x2 = 0; $x2 < $nhat ; $x2++)
                        {
                            $temp_val           = $pr_mat[$p][$x2];
                            $pr_mat[$p][$x2]    = $pr_mat[$p2][$x2];
                            $pr_mat[$p2][$x2]   = $temp_val;
                        }
                    }
                    $fine=1;
                    last;
                }
                
                
            }
            if ($fine==1) {next} ;  # no need to keep looking at this bond for current p and p2 as it is resolved above
            
            ###     STEP 3: ########## IF SELF CANNOT DISTINGUISH PRIORITY  USE R/S RULE ##############
            if ($nors == 1) {next} ;         #printf "not resolved %i %i %i %i %i\n",  $lmat[1][$fpr[$p]],  $lmat[1][$fpr[$p2]], $lmat[1][$i], $lmat[9][$i], $bond_lst_l[1][$i] ;
            $found_4=0;
            if ($bond_lst_l[1][$i] == 4)  {$found_4=1} ;    # four heteroatoms
            if ($bond_lst_l[1][$i] == 3)                    # three heteroatoms
            {
                ########## get attached H
                for ($xi=1 ; $xi <=$numlin ; $xi++)        #   print hydrogens
                {
                    @fields = split (' ', $fil_c[$xi]);
                    if ($xi > $atstart)
                    {
                        if($fields[2] =~ 'H_')
                        {
                            if ($fields[8] == $lmat[1][$i])        # attached carbon
                            {
                                $found_4=1;
                            }
                        }
                    }
                }
                
                ##################
            }
            if ($found_4 == 1)
            {
                #printf "back: %i central: %i  e2: %i e1:%i\n", $lmat[9][$i], $lmat[1][$i], $lmat[1][$fpr[$p]], $lmat[1][$fpr[$p2]];

                for ($xi=$inat ; $xi <=$endat ; $xi++)         #   print hydrogens
                {
                    @fields = split (' ', $fil_c[$xi]);
                    if ($xi > $atstart)
                    {
                        if($fields[0] == $lmat[1][$i])
                        {
                            @at3 = ($fields[2], $fields[3], $fields[4]);
                            if($fields[5] =~ 'C.3')
                            {
                                $found_4=1;
                            }
                        }
                        if($fields[0] == $lmat[9][$i] && $lmat[9][$i] > 0)
                        {
                            #print $lmat[9][$i] $fields0;
                            @at2 = ($fields[2], $fields[3], $fields[4]);
                        }
                        if($fields[0] == $at_n && $lmat[9][$i] == 0)
                        {
                            #print $lmat[9][$i] $fields0;
                            @at2 = ($fields[2], $fields[3], $fields[4]);
                        }
                        
                        if($fields[0] == $lmat[1][$fpr[$p]])
                        {
                            @at1 = ($fields[2], $fields[3], $fields[4]);
                        }
                        if($fields[0] == $lmat[1][$fpr[$p2]])
                        {
                            @at4 = ($fields[2], $fields[3], $fields[4]);
                        }
                        
                    }
                }
                #print "@at1, \n@at2, \n@at3, \n@at4\n";

                getrs(@at1, @at2, @at3, @at4, $det);
                #printf "%f\n", $det;
                #printf "%f\n", $det;
                if ($det < 1)
                {
                    #printf "Clockwise (R / X2). order is: %i %i. \n", $lmat[1][$fpr[$p]],  $lmat[1][$fpr[$p2]];
                    #$fine=1;
                    
                    #last;
                }
                if ($det > 1)
                {
                    #printf "Anti-clockwise (S / X3). order is: %i %i. \n", $lmat[1][$fpr[$p2]],  $lmat[1][$fpr[$p]];
                    #$temp_ch            = $cord[$p];
                    #$cord[$p]           = $cord[$p2];
                    #$cord[$p2]          = $temp_ch;
                    $temp_fp            = $fpr[$p];
                    $fpr[$p]            = $fpr[$p2];
                    $fpr[$p2]           = $temp_fp;
                    
                    for ($x2 = 0; $x2 < $nhat ; $x2++)
                    {
                        $temp_val           = $pr_mat[$p][$x2];
                        $pr_mat[$p][$x2]    = $pr_mat[$p2][$x2];
                        $pr_mat[$p2][$x2]   = $temp_val;
                    }
                    # probably need to set fine=1
                    #$fine=1;
                    #last;
                }
                
            }
            ################### END CIP RULE
            #                                            for ($px = 1 ; $px <= $addchain ; $px++)
            #                                            {
            #                        #                            printf "Atom i: %i Bonds away: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$p]], $x, $pr_mat[$p][$x], $p;
            #                                                    printf "stuff4 A i: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$px]], $cord[$px], $px;
            #                                                    #if ($pr_mat[$p][$x] == 0) {printf "should exit here saying that assignment is arbitrary\n"};
            #                                            }
        }
    }
    #                        for ($px = 1 ; $px <= $addchain ; $px++)
    #                        {
    #    #                            printf "Atom i: %i Bonds away: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$p]], $x, $pr_mat[$p][$x], $p;
    #                                printf "stuff2 A i: %i Priority: %i chain: %i\n", $lmat[1][$fpr[$px]], $cord[$px], $px;
    #                                #if ($pr_mat[$p][$x] == 0) {printf "should exit here saying that assignment is arbitrary\n"};
    #                        }
    
}



