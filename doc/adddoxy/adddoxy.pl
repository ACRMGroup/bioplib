#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    adddoxy
#   File:       adddoxy.pl
#   
#   Version:    V1.0
#   Date:       17.10.14
#   Function:   Create doxygen additional files from comments in C source
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 2014
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#   Reads a block of comments of the form
#   /* Doxygen
#      -------
#      #GROUP    groupname
#      #SUBGROUP subgroupname
#      #FUNCTION functionname()
#      description....
#      description....
#      #KEYFUNCTION functionname()
#      description....
#      description....
#      ...
#      #SUBGROUP subgroupname
#      #FUNCTION functionname()
#      description....
#      description....
#      ...
#      #GROUP    groupname
#      #SUBGROUP subgroupname
#      #FUNCTION functionname()
#      description....
#      description....
#      ...
#   */
#   This is then converted into a set of .dox files for 'extra' 
#   documentation for doxygen.
#
#   A configuration file may be specified with -conf or the default file
#   is adddoxy.conf
#   This is formatted as follows:
#   PACKAGENAME package name   (an overall name for everything)
#   PAGEOFFSET  offset         (the number of .dox files you don't want
#                               to overwrite)
#   OUTDIR      path           (where you want the .dox files to go)
#   GROUP       groupname      (a group for which you want to provide
#                               extra text)
#   GROUPPRE                   (extra text to go before the specified 
#   ...text...                  GROUP)
#   //
#   GROUPPOST                  (extra text to go after the specified 
#   ...text...                  GROUP)
#   //
#   SUBGROUP    subgroupname   (a subgroup for which you want to provide
#                               extra text)
#   SUBGROUPPRE                (extra text to go before the specified 
#   ...text...                  SUBGROUP)
#   //
#   SUBGROUPPOST               (extra text to go after the specified 
#   ...text...                  SUBGROUP)
#   //
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  17.10.14  Original   By: ACRM
#
#*************************************************************************
my $inDoxy = 0;
my $inFunction = 0;
my $tag;
my $details;
my $group = "";
my $subgroup = "";
my $function = "";
my $description = "";
my $keyFunction = 0;

# Take the config file from the command line if specified
# with a -conf= flag - otherwise adddoxy.conf
my $conffile = (defined($::conf)?$::conf:"adddoxy.conf");

$::data        = ();
$::keyData     = ();
$::packageName = "";
$::pageOffset  = 0;
$::outDir      = ".";
ReadConfFile($conffile);        # Read the configuration file

# Read through the C source looking for a comment of the form
# /* Doxygen ...
# */
#
while(<>)
{
    if(/\/\*\s+Doxygen/)        # Found the comment
    {
        $inDoxy    = 1;
        $inFunction = 0;
    }
    elsif($inDoxy)
    {
        chomp;
        s/^\s+//;               # Remove leading spaces
        s/\s+$//;               # Remove trailing spaces

        if(/\*\//)              # End of comment
        {
            ProcessTags($group, $subgroup, $function, $description, $keyFunction);
            $description = "";
            $inDoxy      = 0;
            $inFunction  = 0;
            $keyFunction = 0;
        }
        elsif(/^\#/)            # It's one of our special tags
        {
            # Process existing tag data
            ProcessTags($group, $subgroup, $function, $description, $keyFunction);
            $keyFunction = 0;
            $description = "";

            # Get the tag and the details and remove the # from the tag
            ($tag, $details) = split(/\s+/, $_, 2);
            $tag =~ s/\#//;
            $tag = "\U$tag";
            $inFunction = 0;
            if($tag eq "GROUP") # The #GROUP tag
            {
                $group = $details;
            }
            elsif($tag eq "SUBGROUP") # The #SUBGROUP tag
            {
                $subgroup = $details;
            }
            elsif($tag eq "FUNCTION") # The #FUNCTION tag
            {
                $inFunction  = 1;
                $function    = $details;
                $keyFunction = 0;
            }
            elsif($tag eq "KEYFUNCTION") # The #KEYFUNCTION tag
            {
                $inFunction  = 1;
                $function    = $details;
                $keyFunction = 1;
            }
        }
        elsif($inFunction)
        {
            $description .= " $_";
        }
    }
}

PrintData();


#*************************************************************************
# Process the stored data for the tags for a function
sub ProcessTags
{
    my($group, $subgroup, $function, $description, $keyFunction) = @_;

    if(($description ne "") && 
       ($group       ne "") && 
       ($subgroup    ne "") && 
       ($function    ne ""))
    {
        if(defined($::data{$group}{$subgroup}{$function}) ||
           defined($::keyData{$group}{$subgroup}{$function}))
        {
            print STDERR "$group : $subgroup : $function has been redefined\n";
        }
        
        $description =~ s/^\s+//;
        if($keyFunction)
        {
            $::keyData{$group}{$subgroup}{$function} = $description;
        }
        else
        {
            $::data{$group}{$subgroup}{$function} = $description;
        }
    }
}

#*************************************************************************
# Print any pre list information for a group
sub PrintGroupPre
{
    my($fp, $group) = @_;
    if(defined($::text{$group}{'GROUP'}{'PRE'}))
    {
        print $fp "\n\n";
        print $fp $::text{$group}{'GROUP'}{'PRE'};
        print $fp "\n\n";
    }
}


#*************************************************************************
# Print any post list information for a group
sub PrintGroupPost
{
    my($fp, $group) = @_;
    if(defined($::text{$group}{'GROUP'}{'POST'}))
    {
        print $fp "\n\n";
        print $fp $::text{$group}{'GROUP'}{'POST'};
        print $fp "\n\n";
    }
}

#*************************************************************************
# Print any pre list information for a subgroup
sub PrintSubgroupPre
{
    my($fp, $group, $subgroup) = @_;
    if(defined($::text{$group}{$subgroup}{'PRE'}))
    {
        print $fp "\n\n";
        print $fp $::text{$group}{$subgroup}{'PRE'};
        print $fp "\n\n";
    }
}

#*************************************************************************
# Print any post list information for a group
sub PrintSubgroupPost
{
    my($fp, $group, $subgroup) = @_;
    if(defined($::text{$group}{$subgroup}{'POST'}))
    {
        print $fp "\n\n";
        print $fp $::text{$group}{$subgroup}{'POST'};
        print $fp "\n\n";
    }
}



#*************************************************************************
# Reads the configuration file
sub ReadConfFile
{
    my($conffile) = @_;
    if(open(CONFFILE, $conffile))
    {
        my $group    = "";
        my $subgroup = "";
        my $text     = "";
        while(<CONFFILE>)
        {
            chomp;
            s/^\s+//;
            s/\s+$//;
            my($keyword, $info) = split(/\s+/, $_, 2);
            $keyword = "\U$keyword";
            if($keyword eq "PACKAGENAME")
            {
                $::packageName = $info;
            }
            elsif($keyword eq "PAGEOFFSET")
            {
                $::pageOffset = $info;
            }
            elsif($keyword eq "OUTDIR")
            {
                $::outDir = $info;
            }
            elsif($keyword eq "GROUP")
            {
                $group = $info;
            }
            elsif($keyword eq "SUBGROUP")
            {
                $subgroup = $info;
            }
            elsif($keyword eq "GROUPPRE")
            {
                $::text{$group}{'GROUP'}{'PRE'} = '';
                while(<CONFFILE>)
                {
                    last if(/^\/\//);
                    $::text{$group}{'GROUP'}{'PRE'} .= $_;
                }
            }
            elsif($keyword eq "GROUPPOST")
            {
                $::text{$group}{'GROUP'}{'POST'} = '';
                while(<CONFFILE>)
                {
                    last if(/^\/\//);
                    $::text{$group}{'GROUP'}{'POST'} .= $_;
                }
            }
            elsif($keyword eq "SUBGROUPPRE")
            {
                $::text{$group}{$subgroup}{'PRE'} = '';
                while(<CONFFILE>)
                {
                    last if(/^\/\//);
                    $::text{$group}{$subgroup}{'PRE'} .= $_;
                }
            }
            elsif($keyword eq "SUBGROUPPOST")
            {
                $::text{$group}{$subgroup}{'POST'} = '';
                while(<CONFFILE>)
                {
                    last if(/^\/\//);
                    $::text{$group}{$subgroup}{'POST'} .= $_;
                }
            }
        }
        close CONFFILE;
    }
}


#*************************************************************************
sub PrintData
{
    my $pageCount = $::pageOffset;
    my $pagefile;
    `mkdir -p $::outDir` if(! -d $::outDir);

    foreach my $group (sort keys %::data)
    {
        $pageCount++;
        my $page     = sprintf("page_%02d", $pageCount);
        my $filename = "$::outDir/$page.dox";
        if(open($pagefile, ">$filename"))
        {
            print  $pagefile "/** \n";
            printf $pagefile "\\page page_%02d $group\n", $pageCount;
            printf $pagefile "\\brief\n";
            printf $pagefile "  $::packageName - $group\n\n\n\n\n\n";

            PrintGroupPre($pagefile, $group);

            foreach my $subgroup (sort keys %{$::data{$group}})
            {
                my $title = "$subgroup";
                print $pagefile "\n\n$title\n";
                print $pagefile "-" x length($title);
                print $pagefile "\n\n";

                PrintSubgroupPre($pagefile, $group, $subgroup);

                my $printed = 0;
                foreach my $function (sort keys %{$::keyData{$group}{$subgroup}})
                {
                    if(!$printed)
                    {
                        print $pagefile "\n\n####Key functions\n\n";
                        $printed = 1;
                    }
                    print $pagefile "- $function - $::keyData{$group}{$subgroup}{$function}\n";
                }

                foreach my $function (sort keys %{$::data{$group}{$subgroup}})
                {
                    if($printed)
                    {
                        print $pagefile "\n\n####Other functions\n\n";
                        $printed = 0;
                    }
                    print $pagefile "- $function - $::data{$group}{$subgroup}{$function}\n";
                }

                PrintSubgroupPost($pagefile, $group, $subgroup);
            }

            PrintGroupPost($pagefile, $group);
            print $pagefile "\n\n*\/\n\n";
            close $pagefile;
        }
        else
        {
            print STDERR "Unable to write page file: $filename\n";
            exit 1;
        }
    }
}

