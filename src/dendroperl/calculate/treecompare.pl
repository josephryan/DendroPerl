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

sub symmetric_difference {
    my $tree1 = shift;
    my $tree2 = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in dendropy

    my $ra_t = false_positives_and_negatives($tree1,$tree2,
                                             $is_bipartitions_updated);
    return $ra_t->[0] + $ra_t->[1];
}

sub unweighted_robinson_foulds_distance {
    my $tree1 = shift;
    my $tree2 = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in dendropy

}

sub weighted_robinson_foulds_distance {
}

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

sub false_positives_and_negatives {
    my $reference_tree  = shift;
    my $comparison_tree = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in dendropy

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

sub euclidean_distance {
}

sub find_missing_bipartitions {
    my $reference_tree = shift;
    my $comparison_tree = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in dendropy
}

sub robinson_foulds_distance {
    my $tree1 = shift;
    my $tree2 = shift;
    my $edge_weight_attr = shift;  # default is "length" in dendropy
    $edge_weight_attr = "length" unless ($edge_weight_attr);

}

sub mason_gamer_kellogg_score {
    my $tree1 = shift;
    my $tree2 = shift;
    my $is_bipartitions_updated = shift;  # default is "False" in dendropy
}

sub _get_length_diffs {
}

sub _bipartition_difference {
}

