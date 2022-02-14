#!/usr/bin/perl
#
#Program to rename cyanalib according to IUPAC
#
# ASSUMES ATOM 3 OR ATOM 1 IS ALWAYS THE NITROGEN ATOM
# ASSUMES CARBONYL CARBON IS ALWAYS THE ONE BEFORE CARBONYL OXYGEN UNLESS C-TERM
# ASSUMES THAT CA N and C ARE DEFINED AS BACKBONE ATOMS WITH THOSE NAMES
# 20210923    = M. Mobli (m.mobli@uq.edu.au)

use Getopt::Long;
use File::Basename;
#READ OPTIONS ****************************
GetOptions('i=s' => \$inputfile, 
		   'o=s' => \$outputfile,
			);

#READ LIB FILE **************************
open(LIB, $inputfile) || die "Error: Could not open input library file.\nUsage\n ./fixnom.pl -i input.lib -o output.lib\n ";		# Open the input file
@rlib = <LIB>;										# Read it into an array
close(LIB);
# ***********************
# READ GREEK NAMES
my @gnam  =('A', 'B', 'G', 'D', 'E', 'Z', 'H', 'T', 'I', 'K', 'L', 'M', 'N', 'X', 'O', 'P', 'R', 'S', 'J', 'U', 'F', 'C', 'Y', 'W');
$gind=0;    # index to keep track of names
#********
# READ element types
my @el_t1 = ('PSEUD', 'H_ALI', 'H_AMI', 'H_ARO', 'H_SUL', 'H_OXY', 'C_ALI', 'C_BYL', 'C_ARO', 'C_VIN', 'N_AMI', 'N_AMO', 'O_BYL', 'O_HYD', 'O_EST', 'S_OXY', 'S_RED', 'P_ALI');
my @el_t2 =('Q', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'N', 'N', 'O', 'O', 'O', 'S', 'S', 'P');
my @el_t3 =('0', '1', '1', '1', '1', '1', '6', '6', '6', '6', '7', '7', '8', '8', '8', '16', '16', '15');
my @el_t4 =('0', '1', '1', '1', '1', '1', '6', '6.3', '6.2', '6.1', '7', '7', '8.1', '8', '8.1', '16.1', '16.2', '15');
$el_r=17;


#********

$foundres=$found_n=$found_c=$found_o=0;

	$j=$l=0;
    $i=1;
    $x=0;
$atstart = $hlin = 0;
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
        if ($fields[0] =~ 'ATOMTYPES')
        {
            $numtypes=$fields[1];
            $endattyp = $i+1;
            $stattyp = $i+$numtypes;
        }
        if ($i >= $stattyp && $i <= $endattyp)
        {
            $adat++;
            $atno=$fields[4];
            getattyp($atno,$attyp);
            #printf "Atoms %s **%s %i\n", $fields[1], $attyp, $fields[4];
            $el_t1[17+$adat] = $fields[1];
            $el_t2[17+$adat] = $attyp;
            $el_t3[17+$adat] = $fields[4];
            $el_r++;
        }
		if ($fields[0] =~ 'RESIDUE')
		{

        $foundres = 1;
        $hlin = $i;

        $resnam=$fields[1];

            if ($fields[3]==$fields[5])      # C terminal residue   # Doesn't affect start of frame or end just where CA is
            {
                $cterm=1;
                $at_n=3;
            }
            
            if ($fields[4]==1)      # N terminal residue
            {
                $nterm=1;
                $at_n=1;
                $at_o=$fields[5];
                $at_c = $at_o - 1;

            }
            if($fields[4]==3 && $fields[3]==$fields[5]+1)   # does not work for N or C term residues
            {
                $at_o=$fields[5];
            }
#        $at_o=$fields[5];
        $numang = $fields[2];   # angles
		} # leaving reasidue line

        if ($foundres == 1 && $i > $numang+$hlin)     # after residue statement and after angles.
        {
            if ($atstart == 0 ){ $atstart = $i};        # starts at frist atom after header
                    if ($fields[0] == 3 && $fields[1] == 'N' && $nterm == 0)
                    {
                        $found_n = 1;
                        $at_n=$fields[0] ;
                    }elsif($fields[0] == $at_n && $nterm == 1)    # first nitrogen you hit - shoudl be 1.
                    {
                        $found_n = 1;
                    }
                    if ($fields[0] == $at_o - 1 && $fields[1] == 'C')   # does not work for C-term
                    {
                        $found_c = 1;
                        $at_c = $at_o - 1;
                    }
                    if ($fields[0] == $at_o && $fields[1] == 'O')
                    {
                        $found_o = 1;
                    }
        }
#        printf " %i %i %i %i %i \n", $foundres, $found_n, $found_o, $i, $numang;

        if ($foundres == 1 && $found_n == 1 && $found_o == 0 && $i > $numang)           # doesn't continue past O, cuts out =O and overlap N. - means for C-term the reading frame is the full residue - need to add escape for C=O
        {                                                                               # creating condensed matrix from N amide until Overlap carbonyl oxygen
            unless ($fields[2] =~ 'PSEUD' || $fields[2] =~ 'H_' )
            {
                $l++;
                $temp_crn[$i] = $l;
                $temp_c1[$i] = $fields[0];
                $temp_c2[$i] = $fields[2];
                $temp_c3[$i] = $fields[8];
                $temp_c4[$i] = $fields[9];
                $temp_c5[$i] = $fields[10];
                $temp_c6[$i] = $fields[11];
                $temp_c7[$i] = $fields[12];
                $temp_cx[$i] = $fields[1];
            }
        }
        $fil_c[$i]=$line;           #stores each line in file to array $fil_c
        $i++;                       #counts rows - after row $fields[2] the atoms starts
	}   # end of reading file

$numlin		= $i-1;
$numhead	= $j;
if ($foundres == 0)
{
    
   printf "Error – this is not a residue file. \n";
    exit;

}
#printf "%i %i %i \n", $at_n, $at_ca, $at_c;
#   identify ca - works for all
for ($k=1; $k <= $numlin; $k++)
{
    if ($temp_crn[$k] > 0)
    {
        if($temp_c3[$k] == $at_n ||  $temp_c4[$k] == $at_n || $temp_c5[$k] == $at_n || $temp_c6[$k] == $at_n )  # search condensed array for c_alpha - connected to N and called CA.
        {
            $lnam = length($temp_cx[$k]);
            if ($temp_cx[$k] =~ 'CA' && $lnam == 2)            #ASSUMES THAT CA IS KNOWN AND CONNECTED TO N
            {
                $at_ca = $temp_c1[$k];
                #printf ("%i %s\n", length($temp_cx[$k]), $temp_cx[$k])
            }
            #if($temp_c3[$k] == $at_c ||  $temp_c4[$k] == $at_c || $temp_c5[$k] == $at_c || $temp_c6[$k] == $at_c )
            #{
            #    $at_ca = $temp_c1[$k];
            #}
        }
    }
}

#   find backbone C=O carbon - needed for C-term only
if ($cterm == 1)
{
 for ($k=1; $k <= $numlin; $k++)
 {
    if ($temp_crn[$k] > 0)
    {
        if($temp_c3[$k] == $at_ca ||  $temp_c4[$k] == $at_ca || $temp_c5[$k] == $at_ca || $temp_c6[$k] == $at_ca )  # search condensed array for c_o - connected to CA and called C.
        {
            $lnam = length($temp_cx[$k]);
            if ($temp_cx[$k] =~ 'C' && $lnam == 1)            #ASSUMES THAT CA IS KNOWN AND CONNECTED TO N
            {
                $at_c = $temp_c1[$k];
                #printf ("%i %s\n", length($temp_cx[$k]), $temp_cx[$k])
            }

        }
    }
 }
    for ($k=1; $k <= $numlin; $k++)
    {
       if ($temp_crn[$k] > 0)
       {
           if($temp_c3[$k] == $at_c ||  $temp_c4[$k] == $at_c || $temp_c5[$k] == $at_c || $temp_c6[$k] == $at_c )  # search condensed array for c_o oxygen - connected to C.
           {
               $lnam = length($temp_cx[$k]);
               if ($temp_cx[$k] =~ 'O' && $lnam == 1)            #ASSUMES THAT CA IS KNOWN AND CONNECTED TO N
               {
                   $at_o = $temp_c1[$k];
                   #printf ("%i %s\n", length($temp_cx[$k]), $temp_cx[$k])
               }

           }
       }
    }
}


if ($at_ca < 1 || $at_n < 1 )
{
    
   printf "Error – backbone C alpha and/or N not found. \n";
    exit;

}

#    printf "%i %i %i %i\n", $at_ca, $at_c, $at_n, $at_o;

# CREATE condensed arrays
$nhat = 0;
for ($k=1; $k <= $numlin; $k++)
{
    if ($temp_crn[$k] > 0)
    {
#    printf "b %s %s %s %i %i\n",  $temp_c1[$k],  $temp_c2[$k],  $temp_c3[$k], $i, $k;
        for ($m=1; $m <= $el_r; $m++)
        {
            if ($temp_c2[$k] =~ $el_t1[$m])
            {
#                printf "%s %s %s %s %i t\n", $temp_c2[$k], $el_t1[$m], $el_t2[$m], $el_t3[$m], $m;
                $temp_c7[$k] = $temp_c2[$k] ;
                $temp_c2[$k] = $el_t2[$m];  #delete entry in [column 2] of temp_atom_file and replace with [column 2] of element_types
                #$temp_c8[$k] = $el_t3[$m];  #populate [column 8] of  temp_atom_file with [column 3] of element_types
                next;
            }
        }
        if ($cterm == 1 && $temp_c1[$k] == $at_c){next};    # skip C=O for c-term - this can appear anywhere in the file so needs to be dealt with separately
        if ($cterm == 1 && $temp_c1[$k] == $at_o){next};    # skip C=O for c-term - this can appear anywhere
        $lmat[1][$nhat]=$temp_c1[$k];
        $lmat[2][$nhat]=$temp_c2[$k];
        $lmat[3][$nhat]=$temp_c3[$k];
        $lmat[4][$nhat]=$temp_c4[$k];
        $lmat[5][$nhat]=$temp_c5[$k];
        $lmat[6][$nhat]=$temp_c6[$k];
        $lmat[7][$nhat]=0;  # for later use
        $lmat[8][$nhat]=0;  # for later use
        $lmat[9][$nhat]=0;  # for later use
        $lmat[12][$nhat]=$temp_c7[$k];
#    printf "a %s %s %s %s %s %s %i\n",  $temp_c1[$k],  $temp_c2[$k],  $temp_c3[$k], $temp_c7[$k],  $temp_c5[$k],  $temp_c6[$k],$nhat;
        $nhat++;
    }
}
if ($cterm == 1){$nhat++};  # normally the frame ends at C=O but for C-term the frame ends at the end so need to add an atom to account for the slight shift
    
# condensed array created
$nhat = $nhat-1; # cut out carbonyl carbon start at 1 to cut out backbone N.
# if we know what at_ca is we can skip to here.
#printf "ats %i %i %i %i\n",  $at_n, $at_ca, $at_c, $at_o;
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
    }
    
}

#


# initialise matrix - matrix goes from N amide residue onwards until carbonyl carbon and then stops - carbony carbon must be last residue - can be cut here.

for ($i=1 ; $i < $nhat ; $i++)
{
    $attyp = $lmat[2][$i];          # get atom number  - used for priority
    getatno($attyp, $atno);
    $lmat[10][$i]=$atno;
}
#   FORM CHAINS
#   CREATE PAIRWISE CONNECTIONS GOING BACKWARDS.
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

            for ($j=3; $j <= 6 ; $j++)          # check connected atom for increment $i
            {
                if($lmat[$j][$i] == $ptr)       # something is connected to this point along the chain
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
                        
                }
                
            }
        }
#########END SORTING ROUTINE
}
################## FOR CIP rule need to add a bunch of carbons now. These will need to be checked but not priorities or chained themselves.
### THIS IS NOT WORKING AND Therefore $dum = 0.
$dum=0;
# THIS OVERWRITES THE ENTRY FOR CARBONYL WHICH ISN'T USED ANWYAY
# 8 C_BYL      1.50    0    6   12.0100
# 9 C_ARO      1.60    0    6   12.0100
#10 C_VIN      1.60    0    6   12.0100
#13 O_BYL      1.30   -1    8   16.0000
# Alternatively this would need to be done using information on valency and connected atoms + bond multiplicity
#               Check atom type for each position
#               Determine saturation based on known valence (N=3, C=4, O=2)
#               Determine how many heteroatoms attached
#               Determine if spare valency exists
#                   IF SO
#                       find all atoms connected
#                        for each find if spare valency exists
#               - can't be done as doubel and triple bonds need to be defined at some point
# Currently hard coded for the below
# would be handy for any new atom type if i knew bond multiplicity
# 8 C_BYL      1.50    0    6   12.0100     2   O
# new columns could say number of bonds to - and what atom type number
# could store these in LMAT 14 and 15 for now
# for ($i=1 ; $i < $test ; $i++)  # should be nhat - but not doing this
#
#    if ($lmat[12][$i] =~ 'XXX')
#    {
#        $lmat[1][$nhat+$dum]=1000+$dum;         #atom ID
#        $lmat[2][$nhat+$dum]=$lmat[14][$i];     #atom TYPE
#        $lmat[3][$nhat+$dum]=$lmat[1][$i];      #atom CONNECTION
#        $lmat[4][$nhat+$dum]=0;                 #atom C2
#        $lmat[5][$nhat+$dum]=0;                 #atom C3
#        $lmat[6][$nhat+$dum]=0;                 #atom C4
#        $lmat[7][$nhat+$dum]=$lmat[7][$i];      #atom CHAIN ID
#        $lmat[8][$nhat+$dum]=0;                 #atom C4
#        $lmat[9][$nhat+$dum]=$lmat[1][$i];      #PAIRWISE CONNECTION (BACK TO CA)
#        $lmat[10][$nhat+$dum]=$lmat[15][$i];    #atom NUMBER - FOR PRIORITY
#        $dum++;
#    }
#
#
$test=0;
#if ($test == 1)
#{
#    $bmul_tot=0;        # how many atoms have bond multiplicities
#    $bmul_id[0]=0;      # which atom has bond multiplicity
#    $bmul_n[0]=0;       # how many bond multiplicities (1 for double, 2 for triple etc.)
#    $bmul_an[0][0]=0;   # what is the weight of each of those
#    $bmul_at[0][0]='';  # atom type of each of those
#    # example input file
#    # Tot 2              : How many unique atoms have multiple bonds
#    # ID-MULT 8 2        : Atom 8 has 2 multiple bonds
#    # ID-MULT 12 1       : Atom 12 has 1 multiple bonds
#    # IN 8 2 8 O         : Atom 8 had Double bond to a carbon
#    # IN 8 2 8 O         : Atom 8 also has a Double bond to a carbon
#    # IN 12 3 6 C         : Atom 12 has a triple bond to an oxygen
#    #
#
#    for ($k=1 ; $k < $nhat ; $k++)  # should be nhat - but not doing this
#    {
#        for ($i =1 ; $i < $bmul_tot ; $i++)
#        {
#            if ($lmat[1][$k] == $bmul_id[$i])
#            {
#                for ($j=1 ; $j < $bmul_n[$i] ; $j++)
#                {
#                    $lmat[1][$nhat+$dum]=1000+$dum;         #atom ID
#                    $lmat[2][$nhat+$dum]=$bmul_at[0][0];               #atom TYPE
#                    $lmat[3][$nhat+$dum]=$lmat[1][$i];      #atom CONNECTION
#                    $lmat[4][$nhat+$dum]=0;                 #atom C2
#                    $lmat[5][$nhat+$dum]=0;                 #atom C3
#                    $lmat[6][$nhat+$dum]=0;                 #atom C4
#                    $lmat[7][$nhat+$dum]=$lmat[7][$k];      #atom CHAIN ID
#                    $lmat[8][$nhat+$dum]=0;                 #atom C4
#                    $lmat[9][$nhat+$dum]=$lmat[1][$k];      #PAIRWISE CONNECTION (BACK TO CA)
#                    $lmat[10][$nhat+$dum]=$bmul_an[$i][$j];                #atom NUMBER - FOR PRIORITY
#                    $dum++;
#                }
#            }
#        }
#
#    }
#}


$test=0;
$test=$nhat;

for ($i=1 ; $i < $nhat ; $i++)  # should be nhat - but not doing this
{
    if ($lmat[12][$i] =~ 'C_ARO')
    {
        $lmat[1][$nhat+$dum]=1000+$dum;         #atom ID
        $lmat[2][$nhat+$dum]='C';               #atom TYPE
        $lmat[3][$nhat+$dum]=$lmat[1][$i];      #atom CONNECTION
        $lmat[4][$nhat+$dum]=0;                 #atom C2
        $lmat[5][$nhat+$dum]=0;                 #atom C3
        $lmat[6][$nhat+$dum]=0;                 #atom C4
        $lmat[7][$nhat+$dum]=$lmat[7][$i];      #atom CHAIN ID
        $lmat[8][$nhat+$dum]=0;                 #atom C4
        $lmat[9][$nhat+$dum]=$lmat[1][$i];      #PAIRWISE CONNECTION (BACK TO CA)
        $lmat[10][$nhat+$dum]=6;                #atom NUMBER - FOR PRIORITY
        $dum++;
    }
    if ($lmat[12][$i] =~ 'C_BYL')
    {
        $lmat[1][$nhat+$dum]=1000+$dum;         #atom ID
        $lmat[2][$nhat+$dum]='O';               #atom TYPE
        $lmat[3][$nhat+$dum]=$lmat[1][$i];      #atom CONNECTION
        $lmat[4][$nhat+$dum]=0;                 #atom C2
        $lmat[5][$nhat+$dum]=0;                 #atom C3
        $lmat[6][$nhat+$dum]=0;                 #atom C4
        $lmat[7][$nhat+$dum]=$lmat[7][$i];      #atom CHAIN ID
        $lmat[8][$nhat+$dum]=0;                 #atom C4
        $lmat[9][$nhat+$dum]=$lmat[1][$i];      #PAIRWISE CONNECTION (BACK TO CA)
        $lmat[10][$nhat+$dum]=8;                #atom NUMBER - FOR PRIORITY
        $dum++;
    }
    if ($lmat[12][$i] =~ 'C_VIN')
    {
        $lmat[1][$nhat+$dum]=1000+$dum;         #atom ID
        $lmat[2][$nhat+$dum]='C';               #atom TYPE
        $lmat[3][$nhat+$dum]=$lmat[1][$i];      #atom CONNECTION
        $lmat[4][$nhat+$dum]=0;                 #atom C2
        $lmat[5][$nhat+$dum]=0;                 #atom C3
        $lmat[6][$nhat+$dum]=0;                 #atom C4
        $lmat[7][$nhat+$dum]=$lmat[7][$i];      #atom CHAIN ID
        $lmat[8][$nhat+$dum]=0;                 #atom C4
        $lmat[9][$nhat+$dum]=$lmat[1][$i];      #PAIRWISE CONNECTION (BACK TO CA)
        $lmat[10][$nhat+$dum]=6;                #atom NUMBER - FOR PRIORITY
        $dum++;
    }
    if ($lmat[12][$i] =~ 'O_BYL')
    {
        $lmat[1][$nhat+$dum]=1000+$dum;         #atom ID
        $lmat[2][$nhat+$dum]='C';               #atom TYPE
        $lmat[3][$nhat+$dum]=$lmat[1][$i];      #atom CONNECTION
        $lmat[4][$nhat+$dum]=0;                 #atom C2
        $lmat[5][$nhat+$dum]=0;                 #atom C3
        $lmat[6][$nhat+$dum]=0;                 #atom C4
        $lmat[7][$nhat+$dum]=$lmat[7][$i];      #atom CHAIN ID
        $lmat[8][$nhat+$dum]=0;                 #atom C4
        $lmat[9][$nhat+$dum]=$lmat[1][$i];      #PAIRWISE CONNECTION (BACK TO CA)
        $lmat[10][$nhat+$dum]=6;                #atom NUMBER - FOR PRIORITY
        $dum++;
    }
#    printf "\n";
}
#exit;
$nhat=$nhat+$dum;

#####################################$$$$$$$$$$$$$$$$$$$$$$$$$$
#################  all 1 bond 2 bond 3 bond etc. connectivities for all atoms. - to be used later for assigning priority
for ($i=1 ; $i < $nhat ; $i++)
{
    #printf "%i atoms 1 bond away from %s: \n", $lmat[1][$i], $lmat[12][$i];
    $l=0;
    $master_l[$i][$l] = $lmat[1][$i];   # first entry is itself

    for ($j=1 ; $j < $nhat ; $j++)
    {
        for ($k=3; $k <= 6 ; $k++)          # check connected atom for increment $i
        {
            if ($lmat[1][$i] == $lmat[$k][$j])
            {
     #           printf "%i ", $lmat[1][$j];
                $l++;
                $bond_list[1][$i][$l] = $lmat[1][$j];
                $bond_lst_l[1][$i]=$l;
                $bond_list_ptr[1][$i][$l]=$j;
                $master_l[$i][$l] = $lmat[1][$j];   # add all entries to one bond away
            }
        }
    }
    #printf "\n";
    $master_c[$i]=$l;
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
                        $bond_list[$x+1][$i][$j3]=$bond_list[1][$ptr][$j2];
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
for ($i=1 ; $i < $nhat ; $i++)
{
 # if this is main chain and only entry
    if ($lmat[8][$i] == 1 && $c_ch_at[$lmat[11][$i]] == 1)
    {
        $lmat[13][$i] = sprintf("%s%s ", $lmat[2][$i], $gnam[$lmat[11][$i]]);
        $lmat[13][$i] = sprintf("%-4s", $lmat[13][$i]);
        #if ($lmat[13][$i] =~ 'CA')
        #{
        #    printf "here %i\n", $lmat[9][$i];
        #}
    }else{
        $lmat[13][$i] = sprintf("%s%s%i", $lmat[2][$i], $gnam[$lmat[11][$i]], $lmat[8][$i]);
        $lmat[13][$i] = sprintf("%-4s", $lmat[13][$i]);

    }
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
        if($fields[2] =~ 'H_')
        {
            $l++;
            $hmat[1][$l]=$fields[0];
            $hmat[4][$l]=0;
            $hmat[2][$l]=$fields[8];
            $hmat[7][$l]=$fields[12];       #   pseudo atoms
            $hmat[8][$l]=0;
            $hmat[21][$l]=0;
            if($fields[2] =~ 'H_ARO')       #   FOR QR later
            {
                $hmat[8][$l]=1;
            }
            for ($k=1 ; $k < $nhat ; $k++)
             {
                if($lmat[1][$k] == $fields[8])
                {
                    $hmat[3][$l]=$lmat[11][$k];         # so now I know it is hb or hg etc.
#                    $lmat[15][$k] = $fields[12] ;       #lmat 15 = Q atom associated with those hydrogens
                    $lmat[16][$k]++;                    #lmat 16 = number of H's attached to that carbon total
                    $hmat[5][$l]=$lmat[16][$k];
                }
            }
        if ($fields[8] == 1 && $nterm == 1)
        {
            $allnh++;
            $hmat[21][$l]=$allnh;
            #printf "dol %i\n", $allnh;
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
            if ($lmat[8][$k] == 1 && $c_ch_at[$lmat[11][$k]] == 1)  # stuff done twice beause it is sometimes HB2 HB3 and sometimes HB12 and HB13 - I think
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
                    if ($hmat[5][$i] == $hmat[4][$i])    # if this is the second of a pair
                    {
                         #                            $lmat[30][$k]=$hmat[1][$i];
                         #printf "central (3): %i back(2): %i  h1(4): %i h2(1):%i\n", $lmat[1][$k], $lmat[9][$k], $hmat[1][$i], $lmat[30][$k];
                            # this will be the second hydrogen of the pair being evaluated - the id of the first in the pair is stored in lmat30
                            # put the other in position 1
                        for ($xi=1 ; $xi <=$numlin ; $xi++)         #   print hydrogens
                        {
                            @fields = split (' ', $fil_c[$xi]);
                            if ($xi > $atstart)
                            {
                               # printf "%i %f %f %f\n",  $fields[0], $fields[5], $fields[6], $fields[7];

                                if($fields[0] == $lmat[1][$k])        # PUT SUBSTIUENT WHERE THE CENTRAL ATOM IS AND GIVE LOWEST PRIORITY
                                {
                                    @at3 = ($fields[5], $fields[6], $fields[7]);
                                }
                                if($fields[0] == $at_n && $lmat[9][$k] == 0)    # first nitrogen you hit - shoudl be 1.
                                {
                                    @at2 = ($fields[5], $fields[6], $fields[7]);
                                }
                                
                                if($fields[0] == $lmat[9][$k] && $lmat[9][$k] > 0)
                                {
                                    #print $lmat[9][$i] $fields0;
                                    @at2 = ($fields[5], $fields[6], $fields[7]);
                                }
                                if($fields[0] == $lmat[30][$k])        # first hydrogen atom encountered on that carbon
                                {
                                    @at1 = ($fields[5], $fields[6], $fields[7]);
                                }
                                if($fields[0] == $hmat[1][$i])
                                {
                                    @at4 = ($fields[5], $fields[6], $fields[7]);
                                }
                            }
                        }
                        getrs(@at1, @at2, @at3, @at4, $det);
                        
                        #print "@at1, \n@at2, \n@at3, \n@at4\n $det\n";

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
                    #print "$hmat[5][$i] $hmat[4][$i] $anti\n";

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
                    if ($hmat[5][$i] == $hmat[4][$i])    # if this is the second of a pair
                    {
                        #                            $lmat[30][$k]=$hmat[1][$i];
                        #printf "central (3) : %i back (2): %i  h1(1): %i h2(2):%i\n", $lmat[1][$k], $lmat[9][$k], $hmat[1][$i], $lmat[30][$k];
                        
                        #printf "%i %i  %i %i\n", $lmat[1][$k], $lmat[9][$k], $hmat[1][$i], $lmat[30][$k];
                        for ($xi=1 ; $xi <=$numlin ; $xi++)         #   print hydrogens
                        {
                            @fields = split (' ', $fil_c[$xi]);
                            if ($xi > $atstart)
                            {
                         #       printf "%i %f %f\n",  $fields[0], $fields[6], $fields[7];
                                
                                if($fields[0] == $lmat[1][$k])        # PUT SUBSTIUENT WHERE THE CENTRAL ATOM IS AND GIVE LOWEST PRIORITY
                                {
                                    @at3 = ($fields[5], $fields[6], $fields[7]);
                                }
                                if($fields[0] == $at_n && $lmat[9][$k] == 0)    # first nitrogen you hit - shoudl be 1.
                                {
                                    @at2 = ($fields[5], $fields[6], $fields[7]);
                                }
                                
                                if($fields[0] == $lmat[9][$k] && $lmat[9][$k] > 0)
                                {
                                    #print $lmat[9][$i] $fields0;
                                    @at2 = ($fields[5], $fields[6], $fields[7]);
                                }
                                if($fields[0] == $lmat[30][$k])        # first hydrogen atom encountered on that carbon
                                {
                                    @at1 = ($fields[5], $fields[6], $fields[7]);
                                }
                                if($fields[0] == $hmat[1][$i])
                                {
                                    @at4 = ($fields[5], $fields[6], $fields[7]);
                                }
                            }
                        }
                        #print "@at1, \n@at2, \n@at3, \n@at4\n $det\n";
                        getrs(@at1, @at2, @at3, @at4, $det);
                        if ($det < 1)
                        {
                            #printf "Clockwise. 2\n";
                        }
                        if ($det > 1)       # should only be here on second atom of 2
                        {
                            #printf "Anti-clockwise. 3\n";
                            $anti=1;
                        }
                        
                        
                    }
                    ##################    END   for CIP   ##################
                    #print "$gnam[$lmat[11][$k]], $lmat[8][$k], $hmat[5][$i]+1, \n $anti, $hmat[5][$i], $hmat[4][$i] \n $anti $hmat[5][$i],  $hmat[4][$i]\n";

                    #                    printf "%i %i  \n", $lmat[1][$k], $lmat[9][$k];
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
            $hmat[6][$i] = sprintf("H%i  ", $hmat[21][$i]);     #special case for amide H
        }
        $hmat[20][$i]=1;
        #printf "Yab. %s %i %i\n", $hmat[6][$i], $hmat[1][$i], $allnh;
    }
}
for ($i=1 ; $i <=$hhat ; $i++)        # get pseudoatoms
{
           # printf "%i %s %i\n", $hmat[1][$i], $hmat[2][$i], $hmat[20][$i];
}
##############   CODE FOR PRINTING HYDROGENS END  #################


##############   CODE FOR PRINTING PSEUDO ATOMS  #################

$l=0;
for ($i=1 ; $i <=$numlin ; $i++)        #   print pseudo
{
    @fields = split (' ', $fil_c[$i]);
    if ($i > $atstart)
    {
        if($fields[2] =~ 'PSEUD')
        {
            $l++;
            $qmat[1][$l]=$fields[0];
            $qmat[2][$l]=0;
            $qmat[3][$l]=0;
            $qmat[4][$l]=$fields[12];       # QQ field

             for ($k=1 ; $k <= $hhat ; $k++)
             {
                if($hmat[7][$k] == $fields[0])
                {
                    if ($qmat[2][$l] == 0)              # not first carbon
                    {
                        $lowchain = $hmat[20][$k];
                        $qmat[2][$l]=$hmat[2][$k];         # associated carbon
                    }else{
                        if ($hmat[20][$k] < $lowchain)  # pick carbon with lowest chain
                        {
                            $qmat[2][$l]=$hmat[2][$k];         # associated carbon
                        }
                    }
                    $qmat[20][$l]=$hmat[20][$k];          #same chain, used later for printing
#                    printf "id q %i id c %i chain h %i id h %i k %i\n", $qmat[1][$l], $qmat[2][$l], $qmat[20][$k], $hmat[1][$k], $k ;
#                    printf "Yip %i %i %i.\n", $hmat[1][$k], $qmat[1][$l], $qmat[2][$l];

                    if ($hmat[8][$k] == 1)
                    {
                        $qmat[3][$l]=1;                 # dealing with QR
                        #printf "%i.\n", $qmat[1][$l];

                    }
                }
            }
        }
    }
}
$qhat=$l;

for ($i=1 ; $i <=$qhat ; $i++)        # get pseudoatoms
{
    for ($j=1 ; $j <=$qhat ; $j++)        # get pseudoatoms
    {
        if ($qmat[1][$i] == $qmat[4][$j])
        {
            $qmat[5][$i]=1;
            $qmat[6][$i]=$qmat[2][$j];      #naming carbon
            $qmat[20][$i]=$qmat[20][$j];          #same chain, used later for printing

            if ($qmat[3][$j] == 1)
            {
                $qmat[7][$i]=1;
            }
        }     # this is a QQ or QR
    }
    #printf "Yip %i %i %i.\n", $qmat[1][$i], $qmat[6][$i], $qmat[2][$i];

}

for ($i=1 ; $i <=$qhat ; $i++)        # get pseudoatoms
{
    for ($k=1 ; $k < $nhat ; $k++)
    {


        if ($lmat[1][$k] == $qmat[2][$i] && $qmat[5][$i] == 0)      # not a QQ or QR
        {
            if ($lmat[8][$k] == 1 && $c_ch_at[$lmat[11][$k]] == 1)
            {
                    $qmat[8][$i] = sprintf("Q%s  ", $gnam[$lmat[11][$k]]);
            }elsif ($c_ch_at[$lmat[11][$k]] == 2 && $lmat[16][$k] == 1)      # if 2 chains but only 1 hydrogen on this carbon, then you don't need specifier.
            {
                    $qmat[8][$i] =sprintf("Q%s  ", $gnam[$lmat[11][$k]]);
            }else{
                    $qmat[8][$i] = sprintf("Q%s%i ", $gnam[$lmat[11][$k]], $lmat[8][$k]);
               #printf "Q%s%i \n", $gnam[$lmat[11][$k]], $lmat[8][$k]
            }
        }

        if ($lmat[1][$k] == $qmat[6][$i] && $qmat[5][$i] > 0 && $qmat[7][$i] == 0)      # QQ
        {
            $qmat[8][$i] = sprintf("QQ%s ", $gnam[$lmat[11][$k]]);

        }
        if ($lmat[1][$k] == $qmat[6][$i] && $qmat[5][$i] > 0 && $qmat[7][$i] == 1)      # QR
        {
            #IF FIELDS 2 IS QR THEN JUST OVERWRITE QMAT WITH THE QR stuff.
            $qmat[8][$i] = sprintf("QR  ");

        }
    }
    if (1 == $qmat[2][$i] && $qmat[5][$i] == 0)      # not a QQ or QR
    {
        #printf "nterm %i\n", $qmat[1][$i] ;
        $qmat[8][$i] = sprintf("QHN ");
    }
            #printf "%i %i %i %i %i\n", $qmat[1][$i], $qmat[2][$i], $qmat[20][$i], $qmat[4][$i], $qmat[7][$i];
}
#for ($i=1 ; $i <=$qhat ; $i++)        # get pseudoatoms
#{
#            printf "%i %s %i\n", $qmat[1][$i], $qmat[8][$i], $qmat[20][$i];
#}



#### START PRINTING
    open(OUT, ">$outputfile") || die "Can't open output file";
	my($filename, $directories, $suffix) = fileparse($outputfile);
    $output2 = sprintf("translate_%s", $filename);
    open(OUT2, ">$output2") || die "Can't open output file";
    printf OUT2 "OLD NEW\n";

#        printf OUT "%s",$head[$i];

for ($i=1 ; $i <=$numlin ; $i++)
{
    $pflag=0;
    $hflag=0;
    $qflag=0;

    @fields = split (' ', $fil_c[$i]);
    if ($i > $atstart)
    {

        for ($k=1 ; $k < $nhat ; $k++)
        {
            if ($fields[0] == $lmat[1][$k])
            {
                $pflag=1;
                $pid=$k;
            }
        }
        for ($k=1 ; $k <= $hhat ; $k++)
        {
            if ($fields[0] == $hmat[1][$k])
            {
                $hflag=1;
                $hid=$k;
            }
        }
        for ($k=1 ; $k <= $qhat ; $k++)
        {
            if ($fields[0] == $qmat[1][$k])
            {
                $qflag=1;
                $qid=$k;
                #if ($qmat[7][$i] == 1)  # don't print anything for QR
                #{
                #    $qflag=0;
                #}
            }
        }

        if ($pflag == 1)
        {
            if ($lmat[8][$pid] > 0)
            {
                printf OUT2 "%s %s\n", substr($fil_c[$i],5,4), $lmat[13][$pid];
             if (length($lmat[2][$pid]) == 2)
             {
              substr($fil_c[$i],5,4)=$lmat[13][$pid];
             }else{
              substr($fil_c[$i],5,4)=$lmat[13][$pid];
             }
            }
            printf OUT "%s\n",$fil_c[$i];
                
        }elsif ($hflag == 1)
        {
            if ($hmat[20][$hid] > 0)
            {
                printf OUT2 "%s %s\n", substr($fil_c[$i],5,4), $hmat[6][$hid];
             substr($fil_c[$i],5,4)=$hmat[6][$hid];
            }
            #printf "%s\n %i %i %i\n", $fil_c[$i], $hmat[20][$hid], $hmat[6][$hid], $hid;
            printf OUT "%s\n",$fil_c[$i];
        }elsif ($qflag == 1)
        {
            if ($qmat[20][$qid] > 0 && $qmat[7][$qid] != 1)     # don't change QR
            {
             #printf "%s\n %i %s %i\n", $fil_c[$i], $qmat[20][$qid], $qmat[8][$qid], $qid;
             printf OUT2 "%s %s\n", substr($fil_c[$i],5,4), $qmat[8][$qid];
             substr($fil_c[$i],5,4)=$qmat[8][$qid];
            }
            #printf "%s\n %i %i %i\n", $fil_c[$i], $qmat[20][$qid], $qmat[7][$qid], $qid;

            printf OUT "%s\n",$fil_c[$i];
        }else{
            printf OUT "%s\n",$fil_c[$i];
        }
        
    }else{
        printf OUT "%s\n",$fil_c[$i];
    }
}
#$i++;
#printf OUT "%s",$fil_c[$i];
close (OUT);


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
			#	printf "swap1\n";
                            if ($pr_mat[$p2][0] > $pr_mat[$p][0])
                            {
			#	printf "swap2\n";
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
                        if ($nors == 1) {next} ;
                        #############   START R/S CIP RULE ###################
                            #       printf "not resolved %i %i %i %i %i\n",  $lmat[1][$fpr[$p]],  $lmat[1][$fpr[$p2]], $lmat[1][$i], $lmat[9][$i], $bond_lst_l[1][$i] ;
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

                                for ($xi=1 ; $xi <=$numlin ; $xi++)         #   print hydrogens
                                {
                                    @fields = split (' ', $fil_c[$xi]);
                                    if ($xi > $atstart)
                                    {
                                        #printf "%i %f %f\n",  $fields[0], $fields[6], $fields[7];
                                        
                                        if($fields[0] == $lmat[1][$i])
                                        {
                                            @at3 = ($fields[5], $fields[6], $fields[7]);
                                        }
                                        if($fields[0] == $at_n && $lmat[9][$i] == 0)    # first nitrogen you hit - shoudl be 1. -N is not in loop
                                        {
                                            @at2 = ($fields[5], $fields[6], $fields[7]);
                                        }
                                        
                                        if($fields[0] == $lmat[9][$i] && $lmat[9][$i] > 0)
                                        {
                                            #print $lmat[9][$i] $fields0;
                                            @at2 = ($fields[5], $fields[6], $fields[7]);
                                        }
                                        if($fields[0] == $lmat[1][$fpr[$p]])
                                        {
                                            @at1 = ($fields[5], $fields[6], $fields[7]);
                                        }
                                        if($fields[0] == $lmat[1][$fpr[$p2]])
                                        {
                                            @at4 = ($fields[5], $fields[6], $fields[7]);
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
                        ########### END R/S CIP RULE ##############
                        
                        
                        
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







