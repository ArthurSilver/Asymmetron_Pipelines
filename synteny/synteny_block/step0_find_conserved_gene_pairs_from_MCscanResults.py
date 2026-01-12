import re
import pandas as pd
from collections import defaultdict


def extract_collinear_gene_pairs(file_path, ref_prefix="Alref", target_prefix_pattern=r"B[bfjl]"):
    """
    Extract collinear gene pairs between reference genes and target species genes from a single collinearity file
    :param file_path: Path to collinearity file
    :param ref_prefix: Prefix of reference genome gene ID (e.g., Alref)
    :param target_prefix_pattern: Regular expression for target genome gene ID (e.g., match Bf, Bj, Bl)
    :return: dict, key=reference gene ID, value=target gene ID
    """
    collinear_pairs = dict()
    # Regular expression: match homologous gene pair lines (e.g., "0-  0:	Alref_001268-T1	Bf_018547-T1	 6e-138")
    gene_pair_pattern = re.compile(
        r"\d+-\s*\d+:\s*({0}_\S+)\s+({1}_\S+)".format(ref_prefix, target_prefix_pattern)
    )

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            # Only process lines containing homologous gene pairs (skip parameter, statistics, block title lines)
            match = gene_pair_pattern.search(line)
            if match:
                ref_gene = match.group(1)  # Reference gene ID (e.g., Alref_001268-T1)
                target_gene = match.group(2)  # Target species gene ID (e.g., Bf_018547-T1)
                # Avoid duplicates (if one reference gene corresponds to multiple target genes, take the first reliable pair)
                if ref_gene not in collinear_pairs:
                    collinear_pairs[ref_gene] = target_gene
    return collinear_pairs


def find_conserved_gene_clusters(file_paths):
    """
    Filter reference gene clusters that maintain collinearity across all species
    :param file_paths: List of all collinearity file paths (in target species order)
    :return: dict, key=reference gene ID, value=dict of corresponding gene IDs for each target species
    """
    # Step 1: Extract collinear gene pairs from each file
    all_species_pairs = []
    target_species_names = []  # Record target species names (e.g., Bf, Bj, Bl)
    for file in file_paths:
        # Extract target species name from file name (e.g., extract "Bf" from "Alref.Bf.collinearity.txt")
        species_match = re.search(r"Alref.(B[bfjl]+).collinearity", file)
        if species_match:
            target_species = species_match.group(1)
            target_species_names.append(target_species)
        else:
            raise ValueError(f"File name format error, cannot extract species name: {file}")

        # Extract collinear gene pairs from the file
        gene_pairs = extract_collinear_gene_pairs(file_path=file)
        all_species_pairs.append(gene_pairs)
        print(f"Extracted {len(gene_pairs)} collinear gene pairs from {file}")

    # Step 2: Find intersection of collinear reference genes across all species (core step)
    # First get reference gene set of the first species
    if not all_species_pairs:
        raise ValueError("No collinear gene pairs extracted, please check file paths and format")
    common_ref_genes = set(all_species_pairs[0].keys())

    # Intersect with reference gene sets of subsequent species
    for pairs in all_species_pairs[1:]:
        common_ref_genes.intersection_update(set(pairs.keys()))

    print(f"\nNumber of reference genes maintaining collinearity across all {len(target_species_names)} species: {len(common_ref_genes)}")

    # Step 3: Organize corresponding genes of each collinear reference gene in each target species
    conserved_clusters = defaultdict(dict)
    for ref_gene in common_ref_genes:
        for i, species in enumerate(target_species_names):
            conserved_clusters[ref_gene][species] = all_species_pairs[i][ref_gene]

    return conserved_clusters, target_species_names


def save_conserved_clusters(conserved_clusters, target_species_names, output_file="D:/conserved_collinear_clusters2.csv"):
    """
    Save conserved collinear gene clusters as CSV file
    :param conserved_clusters: Conserved collinear gene cluster dictionary
    :param target_species_names: List of target species names
    :param output_file: Output file path
    """
    # Convert to DataFrame format
    data = []
    for ref_gene, target_genes in conserved_clusters.items():
        row = {"Reference Gene ID": ref_gene}
        # Add corresponding genes for each target species
        for species in target_species_names:
            row[f"{species} Gene ID"] = target_genes[species]
        data.append(row)

    # Save as CSV
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False, encoding="utf-8-sig")
    print(f"\nResults saved to: {output_file}")
    return df


# -------------------------- Main Program Entry --------------------------
if __name__ == "__main__":
    # 1. Configure file paths (modify to your actual file paths, order corresponds to target species)
    file_paths = [
        "filepath/Alref.Bb.collinearity",  # Alref vs Bb
        "filepath/Alref.Bf.collinearity",  # Alref vs Bf
        "filepath/Alref.Bj.collinearity",  # Alref vs Bj
        "filepath/Alref.Bl.collinearity"  # Alref vs Bl
    ]

    # 2. Filter conserved collinear gene clusters
    conserved_clusters, target_species = find_conserved_gene_clusters(file_paths)

    # 3. Save results and print first 10 examples
    result_df = save_conserved_clusters(conserved_clusters, target_species)
    print("\nFirst 10 examples of conserved collinear gene clusters:")
    print(result_df.head(10))