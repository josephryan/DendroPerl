use List::Util qw(sum);
use Math::Complex qw(sqrt);

##############################################################################
# sub symmetric_difference
##############################################################################
# Returns *unweighted* Robinson-Foulds distance between two trees.

# Trees need to share the same |TaxonNamespace| reference. The
# bipartition bitmasks of the trees must be correct for the current tree
# structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
# ``is_bipartitions_updated`` argument must be |False| to force recalculation
# of bipartitions.

# Parameters
# ----------
# $tree1 : |Tree| object
#    The first tree of the two trees being compared. This must share the
#    same |TaxonNamespace| reference as ``tree2`` and must have
#    bipartitions encoded.
# $tree2 : |Tree| object
#    The second tree of the two trees being compared. This must share the
#    same |TaxonNamespace| reference as ``tree1`` and must have
#    bipartitions encoded.
# $is_bipartitions_updated : bool
#    If |False|, then the bipartitions on *both* trees will be updated
#    before comparison. If |True| then the bipartitions will only be
#    calculated for a |Tree| object if they have not been calculated
#    before, either explicitly or implicitly.

# Returns
# ----------
#    The symmetric difference (a.k.a. the unweighted Robinson-Foulds
#    distance) between ``tree1`` and ``tree2``.  (int)
##############################################################################
sub symmetric_difference {
    my $tree1 = shift;
    my $tree2 = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in DendroPy

    my $ra_t = false_positives_and_negatives($tree1,$tree2,
                                             $is_bipartitions_updated);
    return $ra_t->[0] + $ra_t->[1];
}

##############################################################################
# sub unweighted_robinson_foulds_distance
##############################################################################
# Alias for ``symmetric_difference()``.
##############################################################################
sub unweighted_robinson_foulds_distance {
    my $tree1 = shift;
    my $tree2 = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in DendroPy
    return symmetric_difference($tree1,$tree2,$is_bipartitions_updated);
}

##############################################################################
# sub weighted_robinson_foulds_distance
##############################################################################
# Returns *weighted* Robinson-Foulds distance between two trees based on
# ``edge_weight_attr``.
#
# Trees need to share the same |TaxonNamespace| reference. The
# bipartition bitmasks of the trees must be correct for the current tree
# structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
# ``is_bipartitions_updated`` argument must be |False| to force 
# recalculation of bipartitions.

# Parameters
# ----------
# $tree1 : |Tree| object
#    The first tree of the two trees being compared. This must share the
#    same |TaxonNamespace| reference as ``tree2`` and must have
#    bipartitions encoded.
# $tree2 : |Tree| object
#    The second tree of the two trees being compared. This must share the
#    same |TaxonNamespace| reference as ``tree1`` and must have
#    bipartitions encoded.
# $edge_weight_attr : string
#    Name of attribute on edges of trees to be used as the weight.
# $is_bipartitions_updated : bool (1 or 0)
#    If |True|, then the bipartitions on *both* trees will be updated before
#    comparison. If |False| (default) then the bipartitions will only be
#    calculated for a |Tree| object if they have not been calculated
#    before, either explicitly or implicitly.
#
#    Returns
#    -------
#    a float which is 
#    The edge-weighted Robinson-Foulds distance between ``tree1`` and ``tree2``.
##############################################################################
sub weighted_robinson_foulds_distance {
    my $tree1 = shift;
    my $tree2 = shift;
    my $edge_weight_attr = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in DendroPy

    # set default(s)
    $edge_weight_attr = "length" unless ($edge_weight_attr);

    my $df = sub {
        my ($length_diffs) = @_;
        return sum(map { abs($_->[0] - $_->[1]) } @$length_diffs);
    };
    return _bipartition_difference($tree1,$tree2,$df,$edge_weight_attr,
                                   'float',$is_bipartitions_updated);
}

##############################################################################
# sub false_positives_and_negatives
##############################################################################
# Counts and returns number of false positive bipar (bipartitions found in
# ``comparison_tree`` but not in ``reference_tree``) and false negative
# bipartitions (bipartitions found in ``reference_tree`` but not in
# ``comparison_tree``).
# Trees need to share the same |TaxonNamespace| reference. The
# bipartition bitmasks of the trees must be correct for the current tree
# structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
# ``is_bipartitions_updated`` argument must be |False| to force recalculation of
# bipartitions.
# 
# Parameters
# ----------
# $reference_tree : |Tree| object
#     The first tree of the two trees being compared. This must share the
#     same |TaxonNamespace| reference as ``tree2`` and must have
#     bipartitions encoded.
# $comparison_tree : |Tree| object
#     The second tree of the two trees being compared. This must share the
#     same |TaxonNamespace| reference as ``tree1`` and must have
#     bipartitions encoded.
# $is_bipartitions_updated : bool
#     If |True|, then the bipartitions on *both* trees will be updated
#     before comparison. If |False| (default) then the bipartitions
#     will only be calculated for a |Tree| object if they have not been
#     calculated before, either explicitly or implicitly.
# Returns
# -------
# $ra_data
# A pair of integers, with first integer being the number of false
# positives and the second being the number of false negatives.
##############################################################################
sub false_positives_and_negatives {
    my $reference_tree  = shift;
    my $comparison_tree = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in DendroPy

    if ($reference_tree->taxon_namespace ne $comparison_tree->taxon_namespace) {
        die "TaxonNamespaceIdentityError: reference tree namespace ($reference_tree->taxon_namespace) does not equal comparison tree namespace ($comparison_tree->taxon_namespace)";
    }
    if (! $is_bipartitions_updated) {
        $reference_tree->encode_bipartitions();
        $comparison_tree->encode_bipartitions();
    } else {
        if (! defined $reference_tree->bipartition_encoding) {
            $reference_tree->encode_bipartitions();
        }
        if (! defined $comparison_tree->bipartition_encoding) {
            $comparison_tree->encode_bipartitions();
        }
    }
    my %ref_bipartitions = map { $_ => 1 } @{$reference_tree->bipartition_encoding};
    my %comparison_bipartitions = map { $_ => 1 } @{$comparison_tree->bipartition_encoding};
    my @false_positives = grep { ! exists $ref_bipartitions{$_} } keys %comparison_bipartitions;
    my @false_negatives = grep { ! exists $comparison_bipartitions{$_} } keys %ref_bipartitions;
    return [scalar(@false_positives),scalar(@false_negatives)];
}

##############################################################################
# sub euclidean_distance
##############################################################################
#
# Returns the Euclidean distance (a.k.a. Felsenstein's 2004 "branch length
# distance") between two trees based on ``edge_weight_attr``.
#
# Trees need to share the same |TaxonNamespace| reference. The
# bipartition bitmasks of the trees must be correct for the current tree
# structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
# ``is_bipartitions_updated`` argument must be |False| to force recalculation of
# bipartitions.
#
# Parameters
# ----------
# $tree1 : |Tree| object
#        The first tree of the two trees being compared. This must share the
#        same |TaxonNamespace| reference as ``tree2`` and must have
#        bipartitions encoded.
# $tree2 : |Tree| object
#        The second tree of the two trees being compared. This must share the
#        same |TaxonNamespace| reference as ``tree1`` and must have
#        bipartitions encoded.
# $edge_weight_attr : string
#        Name of attribute on edges of trees to be used as the weight.
#    is_bipartitions_updated : bool
#        If |True|, then the bipartitions on *both* trees will be updated
#        before comparison. If |False| (default) then the bipartitions
#        will only be calculated for a |Tree| object if they have not been
#        calculated before, either explicitly or implicitly.
#
# Returns
# -------
# an integer representing
# The Euclidean distance between ``tree1`` and ``tree2``.
##############################################################################
sub euclidean_distance {
    my $tree1 = shift;
    my $tree2 = shift;
    my $edge_weight_attr = shift;
    my $value_type = shift;
    my $is_bipartitions_updated = shift;

    # set default(s)
    $edge_weight_attr = "length" unless ($edge_weight_attr);
    $value_type       = "float"  unless ($value_type);

    my $df = sub {
        my ($length_diffs) = @_;
        return sqrt(sum(map { ($_->[0] - $_->[1]) ** 2 } @$length_diffs));
    };
    return _bipartition_difference($tree1, $tree2, $df,
                                   $edge_weight_attr, $value_type,
                                   $is_bipartitions_updated);
}

##############################################################################
# sub find_missing_bipartitions
##############################################################################
#
# Returns a list of bipartitions that are in ``reference_tree``, but
# not in ``comparison_tree``.
#
# Trees need to share the same |TaxonNamespace| reference. The
# bipartition bitmasks of the trees must be correct for the current tree
# structures (by calling :meth:`Tree.encode_bipartitions()` method) or the
# ``is_bipartitions_updated`` argument must be |False| to force recalculation of
# bipartitions.
#
# Parameters
# ----------
# $reference_tree : |Tree| object
#     The first tree of the two trees being compared. This must share the
#     same |TaxonNamespace| reference as ``tree2`` and must have
#     bipartitions encoded.
# $comparison_tree : |Tree| object
#     The second tree of the two trees being compared. This must share the
#     same |TaxonNamespace| reference as ``tree1`` and must have
#     bipartitions encoded.
# $is_bipartitions_updated : bool
#     If |True|, then the bipartitions on *both* trees will be updated
#     before comparison. If |False| (default) then the bipartitions
#     will only be calculated for a |Tree| object if they have not been
#     calculated before, either explicitly or implicitly.
#
# Returns
# -------
# reference to an array of bipartitions
#     A list of bipartitions that are in the first tree but not in the second.
#
##############################################################################
sub find_missing_bipartitions {
    my $reference_tree = shift;
    my $comparison_tree = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in DendroPy

    my @missing = ();

    if ($reference_tree->taxon_namespace ne $comparison_tree->taxon_namespace) {
        die "TaxonNamespaceIdentityError: reference tree namespace ($reference_tree->taxon_namespace) does not equal comparison tree namespace ($comparison_tree->taxon_namespace)";
    }
    if (!$is_bipartitions_updated) {
        $reference_tree->encode_bipartitions();
        $comparison_tree->encode_bipartitions();
    } else {
        if (not defined $reference_tree->bipartition_encoding()) {
            $reference_tree->encode_bipartitions();
        }
        if (not defined $comparison_tree->bipartition_encoding()) {
            $comparison_tree->encode_bipartitions();
        }
    }
    foreach my $bipartition (@{$reference_tree->bipartition_encoding()}) {
        if (grep { $bipartition eq $_ } @{$comparison_tree->bipartition_encoding}) {
            # The bipartition is present in both trees.
        } else {
            push @missing, $bipartition;
        }
    }
    return \@missing;
}

##############################################################################
# sub robinson_foulds_distance
##############################################################################
#
##############################################################################
sub robinson_foulds_distance {
    my $tree1 = shift;
    my $tree2 = shift;
    my $edge_weight_attr = shift;  # default is "length" in DendroPy

    # set default(s)
    $edge_weight_attr = "length" unless ($edge_weight_attr);

}

##############################################################################
# sub mason_gamer_kellogg_score
##############################################################################
#
##############################################################################
sub mason_gamer_kellogg_score {
    my $tree1 = shift;
    my $tree2 = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in DendroPy
}

##############################################################################
# sub _get_length_diffs
##############################################################################
#  Returns a list of tuples, with the first element of each tuple representing
#  the length of the branch subtending a particular bipartition on ``tree1``, 
#  and the second element the length of the same branch on ``tree2``. If a
#  particular bipartition is found on one tree but not in the other, a value
#  of zero is used for the missing bipartition.
##############################################################################
sub _get_length_diffs {
    my $tree1                   = shift;
    my $tree2                   = shift;
    my $edge_weight_attr        = shift;
    my $value_type              = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in DendroPy
    my $bipartition_length_diff_map = shift; # default is "False" in DendroPy

    # set default(s)
    $edge_weight_attr = "length" unless ($edge_weight_attr);
    $value_type       = 'float'  unless ($value_type);

    my @length_diffs = ();
    my %bipartition_length_diffs = ();
    if ($tree1->taxon_namespace() ne $tree2->taxon_namespace()) {
        die "TaxonNamespaceIdentityError: reference tree namespace ($reference_tree->taxon_namespace) does not equal comparison tree namespace ($comparison_tree->taxon_namespace)";
    }
    if (!$is_bipartitions_updated) {
        $tree1->encode_bipartitions();
        $tree2->encode_bipartitions();
    } else {
        if (!defined($tree1->bipartition_encoding())) {
            $tree1->encode_bipartitions();
        }
        if (!defined($tree2->bipartition_encoding())) {
            $tree2->encode_bipartitions();
        }
    }
    my %tree1_bipartition_edge_map = %{$tree2->bipartition_edge_map()};
    my %tree2_bipartition_edge_map = %{$tree1->bipartition_edge_map()};

    for my $bipartition (keys %tree2_bipartition_edge_map) {
        my $edge = $tree2_bipartition_edge_map{$bipartition};
        my $elen1 = $edge->$edge_weight_attr();
        if (!defined($elen1)) {
            $elen1 = 0;
        }
        my $value1 = $value_type->($elen1);
        my $e2;
        eval { $e2 = delete $tree1_bipartition_edge_map{$bipartition}; };
        my $elen2 = 0;
        if ($@) {
            $elen2 = 0;
        } else {
            $elen2 = $e2->$edge_weight_attr();
            if (!defined($elen2)) {
                if (!defined($e2->tail_node())) {
                    $elen2 = 0.0;
                } else {
                    die "ValueError: Edge length attribute is 'None': Tree: " . $tree2->id() . " ('" . $tree2->label() . "'), Split: " . $bipartition->leafset_as_newick_string($tree2->taxon_namespace());
                }
            }
        }
        my $value2 = $value_type->($elen2);
        push @length_diffs, [$value1, $value2];
        $bipartition_length_diffs{$bipartition} = $length_diffs[-1];
    }
    for my $bipartition (keys %tree1_bipartition_edge_map) {
        my $edge = $tree1_bipartition_edge_map{$bipartition};
        my $elen2 = $edge->$edge_weight_attr();
        if (!defined($elen2)) {
            $elen2 = 0;
        }
        my $value2 = $value_type->($elen2);
        my $e1 = $tree2_bipartition_edge_map{$bipartition};
        my $elen1 = 0.0;
        if (defined($e1)) {
            $elen1 = $e1->$edge_weight_attr();
            if (!defined($elen1)) {
                if (!defined($e1->tail_node())) {
                    $elen1 = 0.0;
                } else {
                    die "ValueError: Edge length attribute is 'None': Tree: $tree1 ('%s'), Split: $bipartition";
                }
            }
        }
        my $value1 = $value_type->($elen1);
        push @length_diffs, [$value1, $value2];
        $bipartition_length_diffs->{$bipartition} = $length_diffs[-1];
    }
#ZZZ LEFT OFF HERE. UNCHECKED chatGPT code below
    if (%bipartition_length_diff_map) {
        return (\@length_diffs, $bipartition_length_diffs);
    } else {
        return \@length_diffs;
    }
}

##############################################################################
# sub _bipartition_difference
##############################################################################
# Returns distance between two trees, each represented by a dictionary of
# bipartitions (as bipartition_mask strings) to edges, using ``dist_fn`` 
# to calculate the distance based on ``edge_weight_attr`` of the edges. 
# ``dist_fn`` is a function that takes a list of pairs of values, where
# the values correspond to the edge lengths of a given bipartition on 
# tree1 and tree2 respectively.
##############################################################################
sub _bipartition_difference {
    my $tree1                   = shift;
    my $tree2                   = shift;
    my $dist_fn                 = shift;
    my $edge_weight_attr        = shift;
    my $value_type              = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in DendroPy

    # set default(s)
    $edge_weight_attr = "length" unless ($edge_weight_attr);
    $value_type       = 'float'  unless ($value_type);

    my $length_diffs  = _get_length_diffs($tree1,$tree2,$edge_weight_attr,
                                          $value_type,$is_bipartitions_updated);
    return dist_fn($length_diffs);
}

##############################################################################
# sub float
##############################################################################
# routine written by Joe Ryan to mimic python's 'float' function
##############################################################################
sub float {
    my $num = shift;

    unless (looks_like_number($num)) {
        die "ValueError: could not convert string to float";
    }

    if (!$num) {
        $num = '0.0';
    } elsif ($num =~ m/\./) {
        $num =~ s/0+$//;
        $num .= 0 if ($num =~ m/\.$/);
    } else {
        $num .= '.0';
    }
    return $num;
}

##############################################################################
# package TreeShapeKernel
##############################################################################
package TreeShapeKernel;

# ZZZ: JOE: I'm not 100% sure what's happening in this chatGPT-generated code
my $TreeShapeKernelNodeCache = __PACKAGE__->namedtuple("_TreeShapeKernelNodeCache",
    ["production", "index", "edge_lengths", "sum_of_square_edge_lengths"]);

sub new {
    my $class = shift;
    my %args  = @_;

    # kernel function
    # sigma=1,
    # gauss_factor=1,
    # decay_factor=0.1,
    my $sigma = delete $args{'sigma'} // 1;
    my $gauss_factor = delete $args{'gauss_factor'} // 1;
    my $decay_factor = delete $args{'decay_factor'} // 0.1;

    my $self = {
        sigma => $sigma,
        gauss_factor => $gauss_factor,
        decay_factor => $decay_factor,
        _tree_cache => {},
    };

    bless $self, $class;

    return $self;
}

# ZZZ: JOE: I'm not 100% sure what's happening in this chatGPT-generated code
sub _TreeShapeKernelNodeCache {
    my ($class, @fields) = @_;
    my $nt = $class;
    my $name = "_TreeShapeKernelNodeCache" . ($nt++);
    my $code = "package $name; use parent 'TreeShapeKernel::_TreeShapeKernelNodeCache'; sub new { my \$class = shift; bless [ \$_[0], \$_[1], \$_[2], \$_[3] ], \$class; } 1;";
    eval $code; die $@ if $@;
    {
        no strict 'refs';
        *{"TreeShapeKernel::$name"} = *{$name . '::'};
    }
    return $name;
}

1;

