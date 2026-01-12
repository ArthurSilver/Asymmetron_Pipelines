import pandas as pd
from collections import defaultdict


def load_tagged_data(input_file):
    """
    Load data with intra-species cluster tags
    :param input_file: Input file path
    :return: Loaded DataFrame (preserves original row order)
    """
    df = pd.read_csv(
        input_file,
        sep='\t',
        dtype=str,  # Avoid scientific notation for position columns
        keep_default_na=False  # Prevent empty values from being automatically replaced with NaN
    )

    # Validate required cluster columns
    required_cluster_cols = [
        'Alref_cluster_tag', 'Bb_cluster_tag',
        'Bf_cluster_tag', 'Bj_cluster_tag', 'Bl_cluster_tag'
    ]
    missing_cols = [col for col in required_cluster_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Input file missing required columns: {', '.join(missing_cols)}")

    # Add original index for each row (used to restore order later)
    df['_original_index'] = range(len(df))
    print(f"Successfully loaded data: {df.shape[0]} rows × {df.shape[1] - 1} columns (including original index column)")
    return df


def get_valid_gene_mask(df):
    """
    Generate "valid gene" mask: exclude genes with cluster = No in any species (only keep genes with valid clusters in all species)
    :param df: DataFrame with intra-species cluster tags
    :return: Boolean mask (True = valid clusters in all species)
    """
    # Masks for valid clusters in each species (exclude No)
    valid_masks = [
        df['Alref_cluster_tag'] != 'No',
        df['Bb_cluster_tag'] != 'No',
        df['Bf_cluster_tag'] != 'No',
        df['Bj_cluster_tag'] != 'No',
        df['Bl_cluster_tag'] != 'No'
    ]

    # Mask for valid in all species (logical AND)
    all_valid_mask = valid_masks[0] & valid_masks[1] & valid_masks[2] & valid_masks[3] & valid_masks[4]
    invalid_count = len(df) - all_valid_mask.sum()

    print(f"\nValid gene filtering:")
    print(f"  Total genes: {len(df)}")
    print(f"  Valid genes (clusters in all species): {all_valid_mask.sum()}")
    print(f"  Invalid genes (No in at least one species): {invalid_count} (will be marked as No)")
    return all_valid_mask


def split_alref_cluster_into_cross_subgroups(df, all_valid_mask):
    """
    Core logic: Split cross-species conserved subgroups within Alref clusters
    Logic: Within the same Alref_cluster, group by (Bb_cluster, Bf_cluster, Bj_cluster, Bl_cluster) quadruple,
          each unique quadruple corresponds to a cross-species major cluster (consistent quadruple = conserved clusters across species)
    :param df: DataFrame with intra-species cluster tags
    :param all_valid_mask: Valid gene mask
    :return: Cross-species major cluster tag dictionary (original index → cross-species major cluster ID)
    """
    # Filter valid gene subset (only analyze genes with clusters in all species)
    valid_df = df[all_valid_mask].copy()

    # Step 1: Group by Alref_cluster_tag, then subdivide by cluster quadruple of other 4 species
    cross_cluster_counter = 0
    cross_tag_map = defaultdict(lambda: 'No')  # Original index → cross-species major cluster ID (default No)
    subgroup_stats = defaultdict(dict)  # Statistics: Alref cluster → quadruple → gene count

    print(f"\nStarting to split Alref clusters and identify cross-species conserved subgroups...")

    # Iterate over each Alref cluster
    for alref_cluster, alref_group in valid_df.groupby('Alref_cluster_tag'):
        # Group by (Bb, Bf, Bj, Bl) cluster quadruple (core splitting logic)
        subgroup_key = ['Bb_cluster_tag', 'Bf_cluster_tag', 'Bj_cluster_tag', 'Bl_cluster_tag']
        subgroup_groups = alref_group.groupby(subgroup_key)

        # Iterate over each subgroup (unique quadruple)
        for (bb_clu, bf_clu, bj_clu, bl_clu), subgroup in subgroup_groups:
            # Assign cross-species major cluster ID to this subgroup
            cross_cluster_counter += 1
            cross_id = f'Cross_{cross_cluster_counter}'

            # Record all genes in this subgroup (map via original index)
            for orig_idx in subgroup['_original_index']:
                cross_tag_map[orig_idx] = cross_id

            # Statistical information (for log output)
            subgroup_size = len(subgroup)
            subgroup_stats[alref_cluster][(bb_clu, bf_clu, bj_clu, bl_clu)] = subgroup_size

            # Print splitting log (for verification)
            print(f"  Alref_{alref_cluster} split subgroup → {cross_id}")
            print(f"    - Number of genes included: {subgroup_size}")
            print(f"    - Corresponding species clusters: Bb={bb_clu}, Bf={bf_clu}, Bj={bj_clu}, Bl={bl_clu}")

    # Output splitting statistics summary
    print(f"\nSplitting completed:")
    print(f"  Original number of valid Alref clusters: {len(valid_df['Alref_cluster_tag'].unique())}")
    print(f"  Number of generated cross-species major clusters: {cross_cluster_counter}")
    print(f"  Details of each Alref cluster split:")
    for alref_clu, subgroups in subgroup_stats.items():
        print(f"    - Alref_{alref_clu}: split into {len(subgroups)} subgroups (total genes: {sum(subgroups.values())})")

    return cross_tag_map, cross_cluster_counter


def add_cross_cluster_column(df, cross_tag_map):
    """
    Add cross-species major cluster column (insert as first column) and restore original row order
    :param df: DataFrame with original index
    :param cross_tag_map: Mapping of original index → cross-species major cluster ID
    :return: Final DataFrame with cross-species major cluster tags
    """
    # Generate cross-species major cluster column (default No, replace valid genes via mapping)
    df['cross_species_cluster_id'] = 'No'
    for orig_idx, cross_id in cross_tag_map.items():
        df.loc[df['_original_index'] == orig_idx, 'cross_species_cluster_id'] = cross_id

    # Restore original row order (sort by _original_index)
    df_final = df.sort_values('_original_index').drop(columns=['_original_index'])

    # Adjust column order: move cross-species major cluster column to first position
    cols = ['cross_species_cluster_id'] + [col for col in df_final.columns if col != 'cross_species_cluster_id']
    df_final = df_final[cols]

    return df_final


def export_result(df_final, output_file='data_with_split_cross_cluster.txt'):
    """
    Export results (tab-separated, preserve original format)
    :param df_final: DataFrame with cross-species major cluster tags
    :param output_file: Output file path
    """
    # Statistics of final results
    cross_stats = df_final['cross_species_cluster_id'].value_counts()
    valid_cross = cross_stats[cross_stats.index.str.startswith('Cross_')]
    no_cross = cross_stats.get('No', 0)

    # Save file
    df_final.to_csv(
        output_file,
        sep='\t',
        index=False,
        encoding='utf-8-sig',
        float_format='%.0f'
    )

    # Print final statistics
    print(f"\n=== Final Result Statistics ===")
    print(f"1. Number of cross-species conserved major clusters: {len(valid_cross)}")
    print(f"2. Total genes in conserved major clusters: {valid_cross.sum()}")
    print(f"3. Number of genes without conserved clusters: {no_cross}")
    print(f"4. Top 5 largest conserved clusters:")
    top5_clusters = valid_cross.head(5)
    for clu, count in top5_clusters.items():
        print(f"   - {clu}: {count} genes")
    print(f"\nResults exported to: {output_file}")
    print(f"Output data dimensions: {df_final.shape[0]} rows × {df_final.shape[1]} columns")


def main():
    # -------------------------- Configuration Parameters (modify according to actual paths) --------------------------
    INPUT_FILE = 'filepath/1.txt'  # Input data (with intra-species cluster tags from step1)
    OUTPUT_FILE = 'filepath/data_with_split_cross_cluster.txt'  # Output path
    # -----------------------------------------------------------------------------------------------------------------

    print("=" * 80)
    print("Cross-species Conserved Major Cluster Tagging Script (supports Alref cluster splitting)")
    print("Core logic: Split subgroups in Alref clusters by (Bb/Bf/Bj/Bl) cluster quadruple, each subgroup is a conserved major cluster")
    print(f"Input file: {INPUT_FILE} | Output file: {OUTPUT_FILE}")
    print("=" * 80)

    try:
        # Step 1: Load data
        print(f"\nStep 1/4: Loading original data...")
        df = load_tagged_data(INPUT_FILE)

        # Step 2: Filter valid genes
        print(f"\nStep 2/4: Filtering valid genes...")
        all_valid_mask = get_valid_gene_mask(df)

        # Step 3: Split Alref clusters and generate cross-species cluster mapping
        print(f"\nStep 3/4: Splitting Alref clusters and identifying conserved subgroups...")
        cross_tag_map, cross_count = split_alref_cluster_into_cross_subgroups(df, all_valid_mask)

        # Step 4: Add cross-species column and export results
        print(f"\nStep 4/4: Generating final data and exporting...")
        df_final = add_cross_cluster_column(df, cross_tag_map)
        export_result(df_final, OUTPUT_FILE)

        print(f"\n=== Processing Completed ===")
        print(f"Total {cross_count} cross-species conserved major clusters generated")

    except Exception as e:
        print(f"\nError during processing: {str(e)}")


if __name__ == "__main__":
    main()