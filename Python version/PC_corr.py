import argparse
import random
import sys
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde, mannwhitneyu, pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score


def C_corr(
    x: np.ndarray, V: np.ndarray, feat_names: list[str], cutoff: float = 0.7
) -> dict:
    """
    Computes the PC-Corr network, filtering features based on PCA loadings and computing correlations.

    Parameters:
    ----------
    x : np.ndarray
        Normalized data matrix of shape (N, M), where N represents the number of samples and M represents the number of features.
    V : np.ndarray
        A 1D array of PCA loadings corresponding to one principal component (length M).
    feat_names : list[str]
        A list containing the feature names, corresponding to the columns of `x`.
    cutoff : float, optional
        The threshold for filtering features based on their absolute PCA loadings (default is 0.7).

    Returns:
    -------
    dict or None
        A dictionary containing the following:
        - 'Edges' : pd.DataFrame
            A DataFrame containing edges of the PC-Corr network, with columns ["Node i", "Node j", "PC-corr(i,j)"].
        - 'Nodes' : pd.DataFrame
            A DataFrame containing nodes (features) and their PCA loadings, with columns ["Node", "Loading (V)"].
        - 'pc_corr' : np.ndarray
            A square matrix of PC-Corr values after applying the cutoff, where each entry represents the PC-Corr value between two features.
        - 'x1' : np.ndarray
            The filtered data matrix (N, M'), containing only the selected features.
        - 'cutoff_f' : float
            The final cutoff used for filtering.

        If no features remain after applying the cutoff, or if no valid edges are found, the function prints a message and returns None.

    Notes:
    ------
    - The function first transforms the PCA loadings (`V`) using a logarithmic transformation and normalization.
    - Features with absolute PCA loadings below the cutoff are removed.
    - Pearson correlation is computed on the remaining features, and the PC-Corr formula is applied.
    - Only edges with absolute PC-Corr values above the cutoff are retained.
    - If no edges remain, the function suggests trying a lower cutoff.
    """

    if len(x.shape) != 2 or len(V.shape) != 1:
        raise ValueError("Input dimensions must be valid: x (NxM), V (M,)")

    V = np.sign(V) * np.log10(1 + np.abs(V) / np.mean(np.abs(V)))
    V = np.sign(V) * (
        (np.abs(V) - np.min(np.abs(V))) / (np.max(np.abs(V)) - np.min(np.abs(V)))
    )

    index = np.abs(V) > cutoff
    if not np.any(index):
        print(f"\nNo features meet the cutoff {cutoff}. Try a lower cutoff.")
        return {}

    removed_features = len(V) - np.sum(index)
    if removed_features > 0:
        print(f"\n{removed_features} feature(s) removed due to |V_feature| < cutoff.")

    x_filtered = x[:, index]
    V_filtered = V[index]
    feat_names_filtered = np.array(feat_names)[index].tolist()

    if x_filtered.shape[1] > 1:
        c = np.corrcoef(x_filtered, rowvar=False)
        pc_corr = np.zeros((len(V_filtered), len(V_filtered)))

        for i in range(len(V_filtered) - 1):
            for j in range(i + 1, len(V_filtered)):
                pc_corr[i, j] = np.sign(c[i, j]) * min(
                    abs(c[i, j]), abs(V_filtered[i]), abs(V_filtered[j])
                )

        pc_corr = pc_corr + pc_corr.T
    else:
        pc_corr = np.zeros((len(V_filtered), len(V_filtered)))

    pc_corr[np.abs(pc_corr) < cutoff] = 0

    mask = np.all(pc_corr == 0, axis=0)
    if np.all(mask):
        print(f"\nNo edges remain after applying cutoff {cutoff}. Try a lower cutoff.")
        return {}

    pc_corr = pc_corr[~mask][:, ~mask]
    sig_PC_corr_Name = np.array(feat_names_filtered)[~mask].tolist()
    x1 = x_filtered[:, ~mask]

    edges_list = []
    for i in range(pc_corr.shape[0] - 1):
        for j in range(i + 1, pc_corr.shape[0]):
            if pc_corr[i, j] != 0:
                edges_list.append(
                    [sig_PC_corr_Name[i], sig_PC_corr_Name[j], pc_corr[i, j]]
                )

    edges = pd.DataFrame(edges_list, columns=["Node i", "Node j", "PC-corr(i,j)"])

    nodes = pd.DataFrame(
        {"Node": sig_PC_corr_Name, "Loading (V)": np.array(V_filtered)[~mask]}
    )

    return {
        "Edges": edges,
        "Nodes": nodes,
        "pc_corr": pc_corr,
        "x1": x1,
        "cutoff_f": cutoff,
    }


def match_V_samplCol(
    colors: list,
    x1: np.ndarray,
    labels: np.ndarray,
    nameLabels: list,
    Nodes: pd.DataFrame,
) -> dict:
    """
    Assigns colors to nodes based on PCA loadings and group-wise mean differences.

    Parameters:
    ----------
    colors : list[str]
        List of colors assigned to different groups.
    x1 : np.ndarray
        Data matrix containing selected features after applying a cutoff.
    labels : np.ndarray
        Array of sample labels corresponding to the rows of x1.
    nameLabels : list[str]
        List of unique group names.
    Nodes : pd.DataFrame
        DataFrame containing feature names and their PCA loadings.

    Returns:
    ----------
    dict
        - "NodeColor": List of assigned colors for each feature (node).
        - "n1_f": Fraction of features assigned to the first color group.
        - "n2_f": Fraction of features assigned to the second color group.

    Notes:
    ------
    - Features are grouped into two main classes (e.g., 'black' and 'red'), with colors assigned based on their
      PCA loadings and mean differences between the two groups.
    - If a feature has a negative loading, it is classified into one of the groups based on its mean difference.
    - The function ensures that features are assigned to the appropriate color groups based on the highest fraction (t-values).
    - Handles ties and cases where one of the groups has no valid entries.
    - Returns a dictionary containing node colors and fractions of assignment.
    """

    n = np.where(np.array(colors) == colors[0])[0][0]  # First color group index
    m = np.where(np.array(colors) == colors[1])[0][0]  # Second color group index

    NodeColor = []
    men1, men2, men3 = [], [], []

    # Compute means of x1 for each feature
    for i in range(x1.shape[1]):
        men1.append(Nodes.iloc[i, 0])  # Feature names
        men2.append(np.mean(x1[labels == nameLabels[n], i]))  # First color group mean
        men3.append(np.mean(x1[labels == nameLabels[m], i]))  # Second color group mean

    men = pd.DataFrame(
        {"Feature": men1, "FirstGroupMean": men2, "SecondGroupMean": men3}
    )

    # Compute differences
    men_dif = np.array(men["FirstGroupMean"]) - np.array(men["SecondGroupMean"])
    m1 = men_dif >= 0
    m2 = men_dif <= 0

    # Identify negative loadings
    b = np.array(Nodes.iloc[:, 1]) <= 0  # Boolean mask for negative loadings
    l1, l2 = np.sum(~b), np.sum(b)

    t1 = np.sum(m1 & ~b) / l1 if l1 > 0 else 0
    t2 = np.sum(m2 & b) / l2 if l2 > 0 else 0
    t3 = np.sum(m1 & b) / l2 if l2 > 0 else 0
    t4 = np.sum(m2 & ~b) / l1 if l1 > 0 else 0

    n1_f, n2_f = np.nan, np.nan

    # Assign colors based on max t values
    if l1 == 0 and t2 >= t3:
        t = t2
        NodeColor = [colors[1] if b[i] else "" for i in range(len(b))]
        n2_f = 1 - t

    elif l1 == 0 and t3 >= t2:
        t = t3
        NodeColor = [colors[0] if b[i] else "" for i in range(len(b))]
        n2_f = 1 - t

    elif l2 == 0 and t1 >= t4:
        t = t1
        NodeColor = [colors[0] if not b[i] else "" for i in range(len(b))]
        n1_f = 1 - t

    elif l2 == 0 and t4 >= t1:
        t = t4
        NodeColor = [colors[1] if not b[i] else "" for i in range(len(b))]
        n1_f = 1 - t

    elif l2 != 0 and l1 != 0 and t1 > t3:
        t = t1
        NodeColor = [colors[0] if not b[i] else colors[1] for i in range(len(b))]
        n1_f, n2_f = 1 - t, 1 - t2

    elif l2 != 0 and l1 != 0 and t3 > t1:
        t = t3
        NodeColor = [colors[1] if not b[i] else colors[0] for i in range(len(b))]
        n1_f, n2_f = 1 - t, 1 - t4

    # Handle tie cases
    if not np.isnan(t1) and not np.isnan(t3) and t1 == 0 and t1 == t3:
        if np.sum(m2 & b) >= np.sum(m2 & ~b):
            NodeColor = [colors[0] if not b[i] else colors[1] for i in range(len(b))]
            n1_f, n2_f = 1 - t1, 1 - t2
        else:
            NodeColor = [colors[1] if not b[i] else colors[0] for i in range(len(b))]
            n1_f, n2_f = 1 - t3, 1 - t4

    if not np.isnan(t2) and not np.isnan(t4) and t2 == 0 and t2 == t4:
        if np.sum(m1 & ~b) >= np.sum(m1 & b):
            NodeColor = [colors[0] if not b[i] else colors[1] for i in range(len(b))]
            n1_f, n2_f = 1 - t1, 1 - t2
        else:
            NodeColor = [colors[1] if not b[i] else colors[0] for i in range(len(b))]
            n1_f, n2_f = 1 - t3, 1 - t4

    return {"NodeColor": NodeColor, "n1_f": n1_f, "n2_f": n2_f}


def get_user_labels():
    while True:
        print(
            "\nIs your data represented by\n[r] ranked labels (labels that are organized according to a progressive order, e.g., different stages of a disease, where Stage 1 < Stage 2 < Stage 3)\n[c] class labels (labels that are not necessarily organized in a progressive order, e.g., Condition A, Condition B, Condition C)? [r/c]:\n"
        )
        u_lab = input("-> ").strip().lower()
        if u_lab in ["r", "c"]:
            break
        print(
            "\nPlease introduce either 'r' for ranked labels or 'c' for class labels\n"
        )

    if u_lab == "r":
        while True:
            print(
                "\nAre the values of your ranked labels\n[d] discrete (Stage 1 < Stage 2 < Stage 3)\n[con] continuous (different times of development of a cell line)? [d/con]:\n"
            )
            u_lab = input("-> ").strip().lower()
            if u_lab in ["d", "con"]:
                break
            print(
                "\nPlease introduce either 'd' for discrete labels or 'con' for continuous labels\n"
            )

    return u_lab


def get_user_normalization():
    while True:
        print(
            "\nThe analysis starts by normalizing or not the dataset.\n\nDo you want to apply:\n[1] no normalization \n[2] a preferred normalization \n[3] automatically all the set of available normalizations? [1/2/3]\n"
        )
        u_norm_opt = input("-> ").strip()
        if u_norm_opt in ["1", "2", "3"]:
            break
        print(
            "\nPlease introduce either 1 for no normalization, 2 for a preferred normalization, or 3 for all the set of available normalizations.\n"
        )

    return int(u_norm_opt)


def get_user_dis_choice():
    while True:
        print("Show sample names in PCA scatter plot? [y/n]\n[y] yes\n[n] no \n")
        dis = input("-> ").strip().lower()

        if dis in ["y", "n"]:
            break
        else:
            print("Please introduce just 'y' or 'n'.")
    return dis


def get_user_cutoff():
    while True:
        try:
            cutoff_input = input(
                "\nSelect a cut-off or a set of cut-offs for generating the PC-corr network [number between 0 and 1]:\n\n"
                "Examples: 0.6 or 0.6, 0.65, 0.7\n\n-> "
            )

            cutoff_values = [float(cutoff) for cutoff in cutoff_input.split(",")]

            if all(0 <= cutoff <= 1 for cutoff in cutoff_values):
                return cutoff_values
            else:
                print(
                    "Please introduce a correct cut-off or a set of cut-offs (in [0,1])."
                )
        except ValueError:
            print("Invalid input. Please enter numbers separated by commas.")


def get_user_aupr(number_el_group: dict, unique_labels: list):
    u_aupr = "r"
    if (
        len(set(number_el_group.values())) == len(number_el_group)
        and len(number_el_group) > 2
    ):
        while True:
            print("\nFor the calculation of AUC and AUPR, choose the positive label:")
            print("[s] Smallest sample group in each comparison")
            print("[l] Largest sample group in each comparison")
            print("[r] Ranked list of possible positive labels")
            u_aupr = input("-> ").strip().lower()
            if u_aupr in ["s", "l", "r"]:
                break
            print("\nPlease choose 's', 'l', or 'r'.")

    u_aupr_r = []
    if u_aupr == "r":
        while True:
            print(
                "\nInput the ranked list of positive labels for AUC and AUPR calculation (comma-separated):"
            )
            shuffled_labels = random.sample(unique_labels, len(unique_labels))
            print("Example:", ", ".join([str(x) for x in shuffled_labels]))
            u_aupr_r = input("-> ").strip().split(",")
            u_aupr_r = [
                label.strip() for label in u_aupr_r if label.strip() in unique_labels
            ]
            if len(u_aupr_r) > 0:
                break
            print("\nPlease enter a valid ranked list of positive labels.")
    return u_aupr, u_aupr_r


def get_user_u_rank(u_lab: str):
    while True:
        if u_lab in ["c", "d"]:
            if u_lab == "c":
                u_rank = input(
                    "\nWould you like to rank the PCA results by \n[p] P-value\n[auc] AUC\n[aupr] AUPR? [p/auc/aupr]:\n\n"
                )
            else:
                u_rank = input(
                    "\nWould you like to rank the PCA results by \n[p] P-value\n[auc] AUC\n[aupr] AUPR\n[pc] Pearson correlation\n[sc] Spearman correlation? [p/auc/aupr/pc/sc]:\n\n"
                )
        else:
            u_rank = input(
                "\nWould you like to rank the PCA results by Pearson or Spearman correlation? [pc/sc]:\n\n"
            )

        if u_rank in ["p", "auc", "aupr", "pc", "sc"]:
            break
        else:
            print("Invalid input. Please enter a valid ranking option.")
    return u_rank


def get_user_u_norm(norms_list: list[str]):
    while True:
        print(
            "\n\nSelect the normalization (for detailed information, see the User guide):\n"
        )

        # Ensure the first sample is not "-"
        poss_norm_case = random.sample(norms_list, len(norms_list))
        while poss_norm_case[0] == "-":
            poss_norm_case = random.sample(norms_list, len(norms_list))

        print(
            f"Examples: {poss_norm_case[0]} or - \nwhere '-' stands for no normalization.\n"
        )
        u_norm_n = input("-> ").strip()

        if u_norm_n in norms_list:
            u_norm = norms_list.index(u_norm_n) + 1
            break
        else:
            print("Please introduce the exact name of the normalization.")
    return u_norm, u_norm_n


def get_user_u_cent():
    while True:
        print(
            "\nCentering version? [y/n]\n[y] yes, centered PCA\n[n] no, non-centered PCA\n"
        )
        u_cent = input("-> ").strip().lower()

        if u_cent in ["y", "n"]:
            break
        else:
            print("Please introduce just 'y' or 'n'.")
    return u_cent


def get_user_u_dim(ncPCA: list):
    while True:
        print("\nSelect the dimension for generating the PC-corr network:\n")

        try:
            u_dim = int(input("-> ").strip())
            if 0 < u_dim <= ncPCA[0].shape[1]:
                u_dim -= 1
                break
            else:
                print("Please introduce an existing dimension.")
        except ValueError:
            print("Invalid input. Please enter a numeric value.")
    return u_dim


def get_user_hide_negative_links():
    while True:
        print(
            "Hide negative links in PC-corr network visualization? [y/n]\n[y] yes\n[n] no \n"
        )
        hide_negative_links = input("-> ").strip().lower()

        if hide_negative_links in ["y", "n"]:
            break
        else:
            print("Please introduce just 'y' or 'n'.")
    return hide_negative_links


def get_user_trustworthiness(u_lab: str):
    if u_lab in ["c", "d"]:
        print(
            "\nAdditionally, you can compute the trustworthiness of p-value, AUC, and AUPR results, "
            "which measures whether the segregation results are discriminative due to randomness "
            "(p-value > 0.05) or if they capture main variability (p-value ≤ 0.05)."
        )

        while True:
            print(
                "Would you like to compute the trustworthiness of the measure of segregation (p-value, AUC, AUPR) but just for the selected case (normalization, dimension, and centering)? [y/n]\n[y] yes\n[n] no \n"
            )
            trustworthiness = input("-> ").strip().lower()

            if trustworthiness in ["y", "n"]:
                break
            else:
                print("Please introduce just 'y' or 'n'.")
        return trustworthiness


def normalize_data(
    x: np.ndarray, u_norm_opt: int
) -> tuple[list[np.ndarray], list[str]]:
    """
    Applies various normalization techniques to the input data matrix.

    Parameters:
    ----------
    x : np.ndarray
        Input data matrix (NxM), where N is the number of samples and M is the number of features.
    u_norm_opt : int
        Option for selecting normalization method:
        - 1: No normalization (returns the original matrix).
        - 2: Interactive selection of a single normalization technique.
        - 3: Applies multiple predefined normalization techniques.

    Returns:
    -------
    norm_matrices : list of np.ndarray
        List containing normalized versions of the input matrix, depending on the selected normalization option.
    norms : list of str
        List of strings indicating the applied normalization techniques.

    Available Normalization Methods:
    --------------------------------
    - "DCS"  : Column-wise division by the sum of each column (Degree-Corrected Scaling).
    - "DRS"  : Row-wise division by the sum of each row (Degree-Regularized Scaling).
    - "LOG"  : Logarithmic transformation (log10(1 + x)).
    - "ZSCORE": Standardization (zero mean, unit variance per column).
    - "PLUS(ABS(MIN))": Shifts all values by adding the absolute minimum.
    - "SQRT" : Element-wise square root transformation.
    - "MANORM": Column-wise division by the mean of each column.

    Notes:
    ------
    - If `u_norm_opt == 1`, only the original matrix is returned with a label "-".
    - If `u_norm_opt == 3`, all normalization techniques are applied sequentially.
    - If `u_norm_opt == 2`, the user is prompted to select a single normalization technique.
    """

    norms = []
    norm_matrices = []

    if u_norm_opt == 1:
        norm_matrices.append(x)
        norms.append("-")
    elif u_norm_opt == 3:
        norm_matrices.append(x / np.sum(x, axis=0, keepdims=True))
        norms.append("DCS")

        norm_matrices.append(x / np.sum(x, axis=1, keepdims=True))
        norms.append("DRS")

        norm_matrices.append(np.log10(1 + x))
        norms.append("LOG")

        norm_matrices.append((x - np.mean(x, axis=0)) / np.std(x, axis=0))
        norms.append("ZSCORE")

        norm_matrices.append(x + abs(np.min(x)))
        norms.append("PLUS(ABS(MIN))")

        norm_matrices.append(np.sqrt(x))
        norms.append("SQRT")

        norm_matrices.append(x / np.mean(x, axis=0, keepdims=True))
        norms.append("MANORM")

    elif u_norm_opt == 2:
        available_norms = [
            "DCS",
            "DRS",
            "LOG",
            "ZSCORE",
            "PLUS(ABS(MIN))",
            "SQRT",
            "MANORM",
        ]
        while True:
            print("\nInput a preferred normalization from the following list:")
            print(", ".join(available_norms))
            u_norm_choice = input("-> ").strip().upper()
            if u_norm_choice in available_norms:
                break
            print("Please introduce the exact name of the normalization.")

        norms.append(u_norm_choice)
        if u_norm_choice == "DCS":
            norm_matrices.append(x / np.sum(x, axis=0, keepdims=True))
        elif u_norm_choice == "DRS":
            norm_matrices.append(x / np.sum(x, axis=1, keepdims=True))
        elif u_norm_choice == "LOG":
            norm_matrices.append(np.log10(1 + x))
        elif u_norm_choice == "ZSCORE":
            norm_matrices.append((x - np.mean(x, axis=0)) / np.std(x, axis=0))
        elif u_norm_choice == "PLUS(ABS(MIN))":
            norm_matrices.append(x + abs(np.min(x)))
        elif u_norm_choice == "SQRT":
            norm_matrices.append(np.sqrt(x))
        elif u_norm_choice == "MANORM":
            norm_matrices.append(x / np.mean(x, axis=0, keepdims=True))

    return norm_matrices, norms


def aupr_evaluation(samp_lab: np.ndarray, scores: np.ndarray, possClass: Any) -> float:
    """
    Computes the Area Under the Precision-Recall Curve (AUPR) for a given set of predictions.

    Parameters:
    ----------
    samp_lab : list
        A list of sample labels, where each label represents the true class of a sample.
    scores : np.ndarray
        A 1D array of predicted scores, where higher values indicate a higher probability of belonging to `possClass`.
    possClass : Any
        The target class for which AUPR is evaluated.

    Returns:
    -------
    float
        The computed AUPR value.
    """
    response = np.array([1 if label == possClass else 0 for label in samp_lab])
    precision, recall, _ = precision_recall_curve(response, scores)
    return auc(recall, precision)


def positive_label_opt(
    u_aupr: str,
    nameLabels1: str,
    nameLabels2: str,
    labels: list,
    u_aupr_r: list[str],
) -> str:
    """
    Determines the positive class label for AUPR (Area Under the Precision-Recall Curve) evaluation.

    Parameters:
    ----------
    u_aupr : str
        Option for determining the positive class:
        - "s": Selects the class with the smaller sample size as the positive class.
        - "l": Selects the class with the larger sample size as the positive class.
        - "r": Selects the class based on its position in the reference ranking list `u_aupr_r`.

    nameLabels1 : str
        Name of the first class label.

    nameLabels2 : str
        Name of the second class label.

    labels : np.ndarray
        Array containing sample labels.

    u_aupr_r : list[str]
        A list defining a ranked order of class labels, used when `u_aupr == "r"`.

    Returns:
    ----------
    str
        The determined positive class label (`possClass`) based on the specified method.

    Notes:
    ------
    -The class appearing earlier in the ranking is considered positive.
    - If a class is not found in `u_aupr_r`, the function defaults to the available class.

    """
    if u_aupr == "s":
        len_n = sum(labels == nameLabels1)
        len_m = sum(labels == nameLabels2)
        possClass = nameLabels1 if len_n < len_m else nameLabels2

    elif u_aupr == "l":
        len_n = sum(labels == nameLabels1)
        len_m = sum(labels == nameLabels2)
        possClass = nameLabels2 if len_n < len_m else nameLabels1

    elif u_aupr == "r":
        idx_n = u_aupr_r.index(nameLabels1) if nameLabels1 in u_aupr_r else None
        idx_m = u_aupr_r.index(nameLabels2) if nameLabels2 in u_aupr_r else None

        if idx_m is None:
            possClass = nameLabels1
        elif idx_n is None:
            possClass = nameLabels2
        else:
            possClass = nameLabels1 if idx_n < idx_m else nameLabels2

    return possClass


def random_permutation_labels(numb_rand: int, labels: list) -> list[list]:
    """
    Generates a specified number of random permutations of the given labels.

    Parameters:
    ----------
    numb_rand : int
        The number of random permutations to generate.

    labels : list
        A list containing the original sample labels.

    Returns:
    ----------
    list[list]
        A list where each element is a randomly shuffled version of `labels`.
    """
    return [random.sample(labels, len(labels)) for _ in range(numb_rand)]


def significance_segregation(
    numb_rand: int,
    PCA: np.ndarray,
    dim: int,
    labels_rp: list[np.ndarray],
    nameLabels: list[str],
    numbLabels: int,
    segr_meas: list[float],
    segr_type: str,
    u_aupr: str,
    u_aupr_r: list[str],
):
    """
    Evaluates the statistical significance of a segregation measure
    (P-value, AUC, or AUPR) through random permutations.

    Parameters:
    ----------
    numb_rand : int
        The number of random label permutations to generate.

    PCA : np.ndarray
        A PCA results matrix where each row corresponds to a sample
        and columns correspond to principal component dimensions.

    dim : int
        The PCA dimension index used for sample segregation.

    labels_rp : list[np.ndarray]
        A list of randomly permuted sample labels for permutation testing.

    nameLabels : list[str]
        A list of unique sample group names.

    numbLabels : int
        The total number of unique sample groups.

    segr_meas : float or list[float]
        The segregation measure(s) computed for the observed data.
        - If `segr_type` is "p", this is a single P-value in a liste.
        - If `segr_type` is "auc" or "aupr", this is a single AUC or AUPR score in a liste.
        - If `segr_type` is "all", this is a list containing [P-value, AUC, AUPR].

    segr_type : str
        The type of segregation measure being evaluated.
        Options:
        - "p"    : P-value (Mann-Whitney U test)
        - "auc"  : Area Under the ROC Curve (AUC)
        - "aupr" : Area Under the Precision-Recall Curve (AUPR)
        - "all"  : Evaluates all three measures (P-value, AUC, and AUPR).

    u_aupr : str
        User-defined option for selecting the positive class in AUPR calculation.

    u_aupr_r : list[str]
        A ranked list of class labels for selecting the positive label.

    Returns:
    ----------
    float or list[float]
        The statistical significance (P-value) of the segregation measure based on
        permutation testing.
        - If `segr_type` is "p", "auc", or "aupr", returns a single float.
        - If `segr_type` is "all", returns a list of three float values:
          [P-value significance, AUC significance, AUPR significance].

    Notes:
    ----------
    - The function generates `numb_rand` permutations of sample labels.
    - It computes the segregation measure for each permutation.
    - The significance is evaluated by comparing the observed measure to
      the distribution of permuted measures.
    - AUC values below 0.5 are adjusted (1 - AUC) for consistency.
    - AUPR scores are computed with flipped scores when AUC is below 0.5.
    """
    rand_segr_meas = np.full((numbLabels * (numbLabels - 1) // 2, numb_rand), np.nan)

    if segr_type == "all":
        rand_segr_meas_mw = np.copy(rand_segr_meas)
        rand_segr_meas_auc = np.copy(rand_segr_meas)
        rand_segr_meas_aupr = np.copy(rand_segr_meas)

    for pr in range(numb_rand):
        labels = labels_rp[pr]
        n, m = 0, 1

        for j in range(numbLabels * (numbLabels - 1) // 2):  # Two-group comparison
            group1 = PCA[np.isin(labels, nameLabels[n]), dim]
            group2 = PCA[np.isin(labels, nameLabels[m]), dim]

            if segr_type == "p":
                rand_segr_meas[j, pr] = mannwhitneyu(
                    group1, group2, alternative="two-sided"
                ).pvalue

            elif segr_type in ["auc", "aupr", "all"]:
                labels_array = np.array(labels)
                samp_lab = np.concatenate(
                    (
                        labels_array[np.isin(labels_array, nameLabels[n])],
                        labels_array[np.isin(labels_array, nameLabels[m])],
                    )
                )
                scores = np.concatenate((group1, group2))

                possClass = positive_label_opt(
                    u_aupr, nameLabels[n], nameLabels[m], labels, u_aupr_r
                )
                negClass = [nameLabels[n], nameLabels[m]]
                negClass.remove(possClass)

                response = np.array(
                    [1 if label == possClass else 0 for label in samp_lab]
                )
                auc_score = roc_auc_score(response, scores)

                if segr_type in ["auc", "all"]:
                    rand_segr_meas_auc[j, pr] = (
                        auc_score if auc_score >= 0.5 else 1 - auc_score
                    )

                if segr_type in ["aupr", "all"]:
                    flip_scores = (
                        2 * np.mean(scores) - scores if auc_score < 0.5 else scores
                    )
                    rand_segr_meas_aupr[j, pr] = aupr_evaluation(
                        samp_lab, flip_scores, possClass
                    )

            m += 1
            if m >= numbLabels:
                n += 1
                m = n + 1

    if segr_type in ["p", "auc", "aupr"]:
        rand_segr_meas_avg = np.mean(rand_segr_meas, axis=1)
    elif segr_type == "all":
        rand_segr_meas_mw_avg = np.mean(rand_segr_meas_mw, axis=1)
        rand_segr_meas_auc_avg = np.mean(rand_segr_meas_auc, axis=1)
        rand_segr_meas_aupr_avg = np.mean(rand_segr_meas_aupr, axis=1)

    if segr_type == "p":
        return (np.sum(rand_segr_meas_avg < segr_meas) + 1) / (numb_rand + 1)
    elif segr_type in ["auc", "aupr"]:
        return (np.sum(rand_segr_meas_avg > segr_meas) + 1) / (numb_rand + 1)
    elif segr_type == "all":
        return [
            (np.sum(rand_segr_meas_mw_avg < segr_meas[0]) + 1) / (numb_rand + 1),
            (np.sum(rand_segr_meas_auc_avg > segr_meas[1]) + 1) / (numb_rand + 1),
            (np.sum(rand_segr_meas_aupr_avg > segr_meas[2]) + 1) / (numb_rand + 1),
        ]


def save_PC_corr_results_to_excel(
    cut_off: list[float], Edges: list[dict], Nodes: list[dict]
) -> None:
    filename = "PC-corr_net_edges-nodes_results.xlsx"

    with pd.ExcelWriter(filename, engine="openpyxl") as writer:
        if isinstance(cut_off, list) and len(cut_off) > 1:
            for edge, node in zip(Edges, Nodes):
                cut_label = "_".join(str(c) for c in cut_off)
                edges_df = edge["data"]
                nodes_df = node["data"]

                if isinstance(edges_df, pd.DataFrame) and isinstance(
                    nodes_df, pd.DataFrame
                ):
                    edges_df.to_excel(
                        writer, sheet_name=f"Edges-{cut_label}", index=False
                    )
                    nodes_df.to_excel(
                        writer, sheet_name=f"Nodes-{cut_label}", index=False
                    )
                else:
                    print(f"Skipping cutoff {cut_label} due to incorrect data format.")
        else:  # Single cutoff case
            cut_label = str(cut_off[0])
            edges_df = Edges[0]["data"]
            nodes_df = Nodes[0]["data"]

            if isinstance(edges_df, pd.DataFrame) and isinstance(
                nodes_df, pd.DataFrame
            ):
                edges_df.to_excel(writer, sheet_name=f"Edges-{cut_label}", index=False)
                nodes_df.to_excel(writer, sheet_name=f"Nodes-{cut_label}", index=False)
            else:
                print("Error: Edges and Nodes are not in DataFrame format!")

    print("\nExcel file saved successfully as 'PC-corr_net_edges-nodes_results.xlsx'.")


def visualize_network(Edges: list[dict], Nodes: list[dict], hide_negative_links: bool):
    """
    Visualizes the PC-Corr network using NetworkX.

    - If two nodes have different colors, their link is drawn as dashed.
    - Positive links are black, negative links are yellow.
    - If hide_negative_links is "y", negative-weight edges are removed.

    Parameters:
    ----------
    Edges : list of dict
        Each dictionary contains a 'cutoff' value and the corresponding edge DataFrame.
    Nodes : list of dict
        Each dictionary contains a 'cutoff' value and the corresponding node DataFrame.
    hide_negative_links : bool
        If True, hides negative-weight edges.
    """

    node_dict = {node_entry["cutoff"]: node_entry["data"] for node_entry in Nodes}

    for edge_entry in Edges:
        cutoff = edge_entry["cutoff"]
        edge_df = edge_entry["data"]

        G = nx.Graph()

        edges_to_draw = []
        edge_colors = []
        edge_styles = []

        for _, row in edge_df.iterrows():
            weight = row["PC-corr(i,j)"]
            if not (hide_negative_links == "y" and weight < 0):
                G.add_edge(row["Node i"], row["Node j"], weight=weight)

        node_df = node_dict.get(cutoff, pd.DataFrame())

        node_colors = {}
        if not node_df.empty:
            for _, row in node_df.iterrows():
                node_colors[row["Node"]] = row["Colour"]

        # Get node positions (force-directed layout)
        pos = nx.spring_layout(G, seed=42, k=10 / G.number_of_nodes())

        node_list = list(G.nodes)
        color_list = [node_colors.get(node, "gray") for node in node_list]
        node_degrees = dict(G.degree())
        node_size_list = [500 + 50 * node_degrees[node] for node in node_list]

        nx.draw_networkx_nodes(
            G, pos, node_color=color_list, node_size=node_size_list, edgecolors="black"
        )

        for node1, node2 in G.edges:
            edges_to_draw.append((node1, node2))

            edge_weight = G[node1][node2]["weight"]
            edge_colors.append("black" if edge_weight >= 0 else "yellow")

            if node_colors.get(node1, "gray") == node_colors.get(node2, "gray"):
                edge_styles.append("solid")
            else:
                edge_styles.append("dashed")

        for edge, color, style in zip(edges_to_draw, edge_colors, edge_styles):
            nx.draw_networkx_edges(
                G,
                pos,
                edgelist=[edge],
                width=1.5,
                alpha=0.6,
                edge_color=color,
                style=style,
            )

        nx.draw_networkx_labels(G, pos, font_size=9, font_color="black")

        line1 = Line2D([0], [0], color="black", lw=2, label="Positive link")
        line2 = Line2D(
            [0], [0], color="black", linestyle="--", lw=2, label="Frustrated link"
        )

        plt.legend(handles=[line1, line2], fontsize=7)

        plt.title(f"c. cutoff={cutoff}", fontsize=18)
        plt.savefig("pc-corr_network.png", dpi=1000)
        plt.show()


def PC_corr_v2(
    df: pd.DataFrame, sample_labels: list[str | int], default_colors=["black", "red"]
) -> tuple[list[dict], list[dict]]:
    """
    Performs PCA-based correlation analysis (PC-Corr) and generates visualizations.

    Parameters:
    ----------
    df : pd.DataFrame
        Input data matrix where:
        - The first row contains feature names.
        - The first column contains sample names.
    sample_labels : list[str | int]
        List of labels corresponding to each sample in df.
    default_colors : list[str], optional
        List of colors for distinguishing sample groups in the PCA plot (default is ['black', 'red']).

    Returns:
    ----------
    tuple
        - Edges (list of dicts): Processed edge data for different cutoffs.
        - Nodes (list of dicts): Processed node data for different cutoffs.

    Visualizes
    ----------
        - PCA scatter plot: Visualization of samples in PCA space.
        - Loadings for each principal component.
        - Significance of loadings.
        - PC-Corr network visualization.
    """

    sample_names = df.iloc[:, 0].tolist()
    feat_names = df.columns[1:].tolist()
    x = df.iloc[:, 1:].values

    labels = sample_labels
    nameLabels: list[Any] = list(set(labels))
    numbLabels = len(nameLabels)

    assert len(sample_labels) == len(
        sample_names
    ), "Error: Length of sample_labels does not match length of sample_names"
    assert not np.isnan(x).any(), "Error: NA values found in the x matrix"
    x = x[:, ~np.all(x == x[0, :], axis=0)]

    u_lab = get_user_labels()
    labl_numb: list[float | int] = []
    if u_lab == "d":
        unique_labels = list(set(sample_labels))
        labl_numb = [unique_labels.index(label) + 1 for label in sample_labels]
    elif u_lab == "con":
        labl_numb = (
            [float(label) for label in sample_labels]
            if isinstance(sample_labels[0], int)
            else []
        )

    u_norm_opt = get_user_normalization()
    norm_matrices, norms = normalize_data(x, u_norm_opt)
    norms_list = norms

    unique_labels = list(set(sample_labels))
    number_el_group = {label: sample_labels.count(label) for label in unique_labels}

    u_aupr, u_aupr_r = get_user_aupr(number_el_group, unique_labels)

    dimension_PCA = min(len(sample_names), len(feat_names))
    pc_nc: list[Any] = [None] * len(norms)
    pc_c: list[Any] = [None] * len(norms)
    ncPCA: list[Any] = [None] * len(norms)
    cPCA: list[Any] = [None] * len(norms)
    explained_nc: list[Any] = [None] * len(norms)
    explained_c: list[Any] = [None] * len(norms)

    mw_ncPCA = np.full(
        (len(norms), len(nameLabels) * (len(nameLabels) - 1) // 2, dimension_PCA),
        np.nan,
    )
    mw_cPCA = np.copy(mw_ncPCA)
    AUC_nc = np.copy(mw_ncPCA)
    AUC_c = np.copy(AUC_nc)
    AUPR_nc = np.copy(AUC_nc)
    AUPR_c = np.copy(AUC_nc)
    rank_pears_corr_ncPCA = np.full((len(norms), dimension_PCA, 1), np.nan)
    rank_pears_corr_cPCA = np.copy(rank_pears_corr_ncPCA)
    rank_spear_corr_ncPCA = np.copy(rank_pears_corr_ncPCA)
    rank_spear_corr_cPCA = np.copy(rank_pears_corr_ncPCA)

    for i in range(len(norms)):
        # Non-centered PCA
        pca_nc = PCA(svd_solver="full", whiten=False)
        ncPCA[i] = pca_nc.fit_transform(norm_matrices[i].copy())
        pc_nc[i] = pca_nc.components_
        explained_nc[i] = 100 * pca_nc.explained_variance_ratio_

        # Centered PCA
        normCenter = norm_matrices[i] - np.mean(norm_matrices[i], axis=0)[np.newaxis, :]

        pca_c = PCA(svd_solver="full", whiten=False)
        cPCA[i] = pca_c.fit_transform(normCenter)
        pc_c[i] = pca_c.components_
        explained_c[i] = 100 * pca_c.explained_variance_ratio_

        for k in range(ncPCA[i].shape[1]):
            if u_lab in ["c", "d"]:
                n, m = 0, 1
                for j in range(
                    len(nameLabels) * (len(nameLabels) - 1) // 2
                ):  # Two-group comparison
                    # Mann-Whitney test
                    group1_nc = ncPCA[i][np.isin(labels, nameLabels[n]), k]
                    group2_nc = ncPCA[i][np.isin(labels, nameLabels[m]), k]
                    mw_ncPCA[i, j, k] = mannwhitneyu(
                        group1_nc, group2_nc, alternative="two-sided"
                    ).pvalue

                    group1_c = cPCA[i][np.isin(labels, nameLabels[n]), k]
                    group2_c = cPCA[i][np.isin(labels, nameLabels[m]), k]
                    mw_cPCA[i, j, k] = mannwhitneyu(
                        group1_c, group2_c, alternative="two-sided"
                    ).pvalue

                    labels = np.array(labels)
                    samp_lab = np.concatenate(
                        (
                            labels[np.isin(labels, nameLabels[n])],
                            labels[np.isin(labels, nameLabels[m])],
                        )
                    )
                    scores_nc = np.concatenate((group1_nc, group2_nc))
                    scores_c = np.concatenate((group1_c, group2_c))

                    possClass = positive_label_opt(
                        u_aupr, nameLabels[n], nameLabels[m], labels, u_aupr_r
                    )
                    negClass = [nameLabels[n], nameLabels[m]]
                    negClass.remove(possClass)

                    response = np.array(
                        [1 if label == possClass else 0 for label in samp_lab]
                    )

                    # Compute AUC
                    AUC_nc[i, j, k] = roc_auc_score(response, scores_nc)
                    AUC_c[i, j, k] = roc_auc_score(response, scores_c)

                    if AUC_nc[i, j, k] < 0.5:
                        AUC_nc[i, j, k] = 1 - AUC_nc[i, j, k]
                        flip_scores_nc = 2 * np.mean(scores_nc) - scores_nc
                        AUPR_nc[i, j, k] = aupr_evaluation(
                            samp_lab, flip_scores_nc, possClass
                        )
                    else:
                        AUPR_nc[i, j, k] = aupr_evaluation(
                            samp_lab, scores_nc, possClass
                        )

                    if AUC_c[i, j, k] < 0.5:
                        AUC_c[i, j, k] = 1 - AUC_c[i, j, k]
                        flip_scores_c = 2 * np.mean(scores_c) - scores_c
                        AUPR_c[i, j, k] = aupr_evaluation(
                            samp_lab, flip_scores_c, possClass
                        )
                    else:
                        AUPR_c[i, j, k] = aupr_evaluation(samp_lab, scores_c, possClass)

                    m += 1
                    if m >= len(nameLabels):
                        n += 1
                        m = n + 1

                if u_lab == "d":
                    rank_pears_corr_ncPCA[i, k, 0] = pearsonr(
                        ncPCA[i][:, k], labl_numb
                    )[0]
                    rank_pears_corr_cPCA[i, k, 0] = pearsonr(cPCA[i][:, k], labl_numb)[
                        0
                    ]
                    rank_spear_corr_ncPCA[i, k, 0] = spearmanr(
                        ncPCA[i][:, k], labl_numb
                    )[0]
                    rank_spear_corr_cPCA[i, k, 0] = spearmanr(cPCA[i][:, k], labl_numb)[
                        0
                    ]
            else:
                rank_pears_corr_ncPCA[i, k, 0] = pearsonr(ncPCA[i][:, k], labl_numb)[0]
                rank_pears_corr_cPCA[i, k, 0] = pearsonr(cPCA[i][:, k], labl_numb)[0]
                rank_spear_corr_ncPCA[i, k, 0] = spearmanr(ncPCA[i][:, k], labl_numb)[0]
                rank_spear_corr_cPCA[i, k, 0] = spearmanr(cPCA[i][:, k], labl_numb)[0]

    dim = np.arange(1, len(ncPCA[0]) + 1)
    dim = np.repeat(dim, len(norms))
    dim = np.concatenate((dim, dim))
    centred = np.repeat("yes", len(ncPCA[0]) * len(norms))
    non_centred = np.repeat("no", len(ncPCA[0]) * len(norms))
    centr = np.concatenate((non_centred, centred))

    norms = np.tile(norms, len(ncPCA[0]) * 2)

    # Convert explained variance lists to NumPy arrays
    explained_c_array = np.concatenate([np.array(ec) for ec in explained_c])
    explained_nc_array = np.concatenate([np.array(enc) for enc in explained_nc])

    # Ensure reshaping works by checking array sizes
    variance_c = explained_c_array.flatten()
    variance_nc = explained_nc_array.flatten()

    variance = np.concatenate((variance_nc, variance_c))

    # Ensure all variables are NumPy arrays for easy calculations
    norms = np.array(norms)
    centr = np.array(centr)
    dim = np.array(dim)
    variance = np.array(variance)

    results = pd.DataFrame(
        {
            "P-value": np.concatenate(
                [np.array(mw_ncPCA).flatten(), np.array(mw_cPCA).flatten()]
            ),
            "AUC": np.concatenate(
                [np.array(AUC_nc).flatten(), np.array(AUC_c).flatten()]
            ),
            "AUPR": np.concatenate(
                [np.array(AUPR_nc).flatten(), np.array(AUPR_c).flatten()]
            ),
            "pears": np.concatenate(
                [
                    np.array(rank_pears_corr_ncPCA).flatten(),
                    np.array(rank_pears_corr_cPCA).flatten(),
                ]
            ),
            "spear": np.concatenate(
                [
                    np.array(rank_spear_corr_ncPCA).flatten(),
                    np.array(rank_spear_corr_cPCA).flatten(),
                ]
            ),
            "Norm": np.array(norms).flatten(),
            "Centering": np.array(centr).flatten(),
            "Dim": np.array(dim).flatten(),
            "expl Var": np.array(variance).flatten(),
        }
    )

    print(
        "\nThe best discrimination in a PCA result (combination of normalization and centering) and along one dimension (principal component, PCn) can be assessed by different evaluators.\n"
    )

    u_rank = get_user_u_rank(u_lab)

    print(f"{u_rank=}")

    ranking_column = {
        "p": "P-value",
        "auc": "AUC",
        "aupr": "AUPR",
        "pc": "pears",
        "sc": "spear",
    }[u_rank]

    ascending = True if u_rank == "p" else False
    results = results.sort_values(by=ranking_column, ascending=ascending)

    # Display a subset of the best results
    if u_rank in ["p", "auc", "aupr"]:
        threshold = 0.05 if u_rank == "p" else 0.7
        sub_results = (
            results[results[ranking_column] < threshold]
            if u_rank == "p"
            else results[results[ranking_column] >= threshold]
        )
    else:
        sub_results = results[abs(results[ranking_column]) >= 0.6]

    if sub_results.empty:
        print(
            f"\nThe table below shows the best results ({ranking_column} >= {0.6 if u_rank in ['auc', 'aupr'] else 0.05}), ranked from most to least discriminative.\n"
        )
        print(results.head())  # Show top few if none meet threshold
    else:
        print(
            f"\nThe table below shows the best results ({ranking_column} >= {0.6 if u_rank in ['auc', 'aupr'] else 0.05}), ranked from most to least discriminative.\n"
        )
        print(sub_results)

    filename = "results.xlsx"
    results.to_excel(filename, sheet_name="PCA results", index=False)
    print(f"\nAll results saved to {filename}, ranked by {ranking_column}.\n")

    print(
        "\nAfter seeing the ranked results, you can choose any combination of normalization, centering, dimension (PCn), and cut-off for the network to obtain the PC-corr network according to your interest (or need).\n"
    )

    if u_norm_opt in [1, 2]:
        u_norm = 1
        u_norm_n = norms_list[0]
    else:
        u_norm, u_norm_n = get_user_u_norm(norms_list)

    u_cent = get_user_u_cent()

    u_dim = get_user_u_dim(ncPCA)

    if u_cent == "y":
        PCAA = cPCA
        if u_lab == "c":
            if u_rank == "p":
                mw_PCA = mw_cPCA
                evaluat = "P-value"
            elif u_rank == "auc":
                mw_PCA = AUC_c
                evaluat = "AUC"
            elif u_rank == "aupr":
                mw_PCA = AUPR_c
                evaluat = "AUPR"

        elif u_lab == "d":
            mw_PCA3 = np.transpose(rank_pears_corr_cPCA, (0, 2, 1))
            mw_PCA4 = np.transpose(rank_spear_corr_cPCA, (0, 2, 1))

            if u_rank == "p":
                mw_PCA = mw_cPCA
                evaluat = "P-value"
            elif u_rank == "auc":
                mw_PCA = AUC_c
                evaluat = "AUC"
            elif u_rank == "aupr":
                mw_PCA = AUPR_c
                evaluat = "AUPR"
            elif u_rank == "pc":
                mw_PCA = mw_PCA3
                mw_PCA3 = AUPR_c
                evaluat = "Pearson correlation"
            elif u_rank == "sc":
                mw_PCA = mw_PCA4
                mw_PCA3 = AUPR_c
                evaluat = "Spearman correlation"

        else:
            if u_rank == "pc":
                mw_PCA = np.transpose(rank_pears_corr_cPCA, (0, 2, 1))
                evaluat = "Pearson correlation"
            else:
                mw_PCA = np.transpose(rank_spear_corr_cPCA, (0, 2, 1))
                evaluat = "Spearman correlation"

        pc = pc_c
        ttl = "centred PCA"
        explained = explained_c

    elif u_cent == "n":
        PCAA = ncPCA
        if u_lab == "c":
            if u_rank == "p":
                mw_PCA = mw_ncPCA
                evaluat = "P-value"
            elif u_rank == "auc":
                mw_PCA = AUC_nc
                evaluat = "AUC"
            elif u_rank == "aupr":
                mw_PCA = AUPR_nc
                evaluat = "AUPR"

        elif u_lab == "d":
            mw_PCA3 = np.transpose(rank_pears_corr_ncPCA, (0, 2, 1))
            mw_PCA4 = np.transpose(rank_spear_corr_ncPCA, (0, 2, 1))

            if u_rank == "p":
                mw_PCA = mw_ncPCA
                evaluat = "P-value"
            elif u_rank == "auc":
                mw_PCA = AUC_nc
                evaluat = "AUC"
            elif u_rank == "aupr":
                mw_PCA = AUPR_nc
                evaluat = "AUPR"
            elif u_rank == "pc":
                mw_PCA = mw_PCA3
                mw_PCA3 = AUPR_nc
                evaluat = "Pearson correlation"
            elif u_rank == "sc":
                mw_PCA = mw_PCA4
                mw_PCA3 = AUPR_nc
                evaluat = "Spearman correlation"

        else:
            if u_rank == "pc":
                mw_PCA = np.transpose(rank_pears_corr_ncPCA, (0, 2, 1))
                evaluat = "Pearson correlation"
            else:
                mw_PCA = np.transpose(rank_spear_corr_ncPCA, (0, 2, 1))
                evaluat = "Spearman correlation"

        pc = pc_nc
        ttl = "non-centred PCA"
        explained = explained_nc

    ind2 = None

    if u_lab in ["c", "d"]:
        if u_rank == "p":

            u_norm = min(u_norm, mw_PCA.shape[0] - 1)

            val = np.min(mw_PCA[u_norm, 0, :])  # Indexing in Python is 0-based
            ind1 = np.argmin(mw_PCA[u_norm, 0, :])  # Get index of min value

            if ind1 == u_dim:
                if ind1 == mw_PCA.shape[2] - 1:  # Last index case
                    val = np.min(mw_PCA[u_norm, 0, :ind1])  # Exclude last index
                    ind2 = np.argmin(mw_PCA[u_norm, 0, :ind1])
                elif ind1 == 0:  # First index case
                    val = np.min(mw_PCA[u_norm, 0, 1:])  # Exclude first index
                    ind2 = np.argmin(mw_PCA[u_norm, 0, 1:]) + 1  # Adjust for offset
                else:  # General case (middle indices)
                    mask = np.arange(mw_PCA.shape[2]) != ind1  # Mask to exclude ind1
                    val = np.min(mw_PCA[u_norm, 0, mask])
                    ind2 = np.argmin(mw_PCA[u_norm, 0, mask])

        elif u_rank in ["pc", "sc"]:  # Pearson or Spearman correlation

            val = np.max(np.abs(mw_PCA[u_norm, 0, :]))
            ind1 = np.argmax(np.abs(mw_PCA[u_norm, 0, :]))

            if ind1 == u_dim:
                if ind1 == mw_PCA.shape[2] - 1:  # Last index case
                    val = np.max(np.abs(mw_PCA[u_norm, 0, :ind1]))
                    ind2 = np.argmax(np.abs(mw_PCA[u_norm, 0, :ind1]))
                elif ind1 == 0:  # First index case
                    val = np.max(np.abs(mw_PCA[u_norm, 0, 1:]))
                    ind2 = (
                        np.argmax(np.abs(mw_PCA[u_norm, 0, 1:])) + 1
                    )  # Adjust for offset
                else:  # General case (middle indices)
                    mask = np.arange(mw_PCA.shape[2]) != ind1  # Mask to exclude ind1
                    val = np.max(np.abs(mw_PCA[u_norm, 0, mask]))
                    ind2 = np.argmax(np.abs(mw_PCA[u_norm, 0, mask]))

        elif u_rank in ["auc", "aupr"]:  # AUC and AUPR

            val = np.max(mw_PCA[u_norm, 0, :])
            ind1 = np.argmax(mw_PCA[u_norm, 0, :])

            if ind1 == u_dim:
                if ind1 == mw_PCA.shape[2] - 1:
                    val = np.max(mw_PCA[u_norm, 0, :ind1])
                    ind2 = np.argmax(mw_PCA[u_norm, 0, :ind1])
                elif ind1 == 0:
                    val = np.max(mw_PCA[u_norm, 0, 1:])
                    ind2 = np.argmax(mw_PCA[u_norm, 0, 1:]) + 1
                else:
                    mask = np.arange(mw_PCA.shape[2]) != ind1
                    val = np.max(mw_PCA[u_norm, 0, mask])
                    ind2 = np.argmax(mw_PCA[u_norm, 0, mask])

    else:
        val = np.max(np.abs(mw_PCA[u_norm, 0, :]))
        ind1 = np.argmax(np.abs(mw_PCA[u_norm, 0, :]))

        if ind1 == u_dim:
            if ind1 == mw_PCA.shape[2] - 1:
                val = np.max(np.abs(mw_PCA[u_norm, 0, :ind1]))
                ind2 = np.argmax(np.abs(mw_PCA[u_norm, 0, :ind1]))
            elif ind1 == 0:
                val = np.max(np.abs(mw_PCA[u_norm, 0, 1:]))
                ind2 = np.argmax(np.abs(mw_PCA[u_norm, 0, 1:])) + 1
            else:
                mask = np.arange(mw_PCA.shape[2]) != ind1
                val = np.max(np.abs(mw_PCA[u_norm, 0, mask]))
                ind2 = np.argmax(np.abs(mw_PCA[u_norm, 0, mask]))

    if ind2 is None:
        ind = ind1
    else:
        ind = ind2 + 1 if ind2 >= ind1 else ind2

    # Determine dimensions for scatter plot
    if u_rank == "p":
        if val <= mw_PCA[u_norm, 0, u_dim]:
            dim1, dim2 = ind, u_dim
        else:
            dim1, dim2 = u_dim, ind
    else:
        if abs(val) >= abs(mw_PCA[u_norm, 0, u_dim]):
            dim1, dim2 = ind, u_dim
        else:
            dim1, dim2 = u_dim, ind

    # Assign colors to groups
    group_medians = [
        np.median(PCAA[u_norm][labels == name, dim1]) for name in nameLabels
    ]
    sorted_indices = np.argsort(group_medians)
    colors = [default_colors[0]] * numbLabels
    colors[sorted_indices[-1]] = default_colors[1]

    dis = get_user_dis_choice()

    # Create density plots for the x-axis
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, name in enumerate(nameLabels):
        grp_data = PCAA[u_norm][labels == name, dim1]
        density = gaussian_kde(grp_data)
        x_vals = np.linspace(grp_data.min(), grp_data.max(), 100)
        y_vals = density(x_vals)
        ax.plot(x_vals, y_vals, color=colors[i])

    ax.set_xlabel(f"PC{dim1+1}")
    ax.set_ylabel("Density")
    plt.savefig("pca_density_x.png", dpi=300)

    # Create density plots for the y-axis
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, name in enumerate(nameLabels):
        grp_data = PCAA[u_norm][labels == name, dim2]
        density = gaussian_kde(grp_data)
        x_vals = np.linspace(grp_data.min(), grp_data.max(), 100)
        y_vals = density(x_vals)
        ax.plot(x_vals, y_vals, color=colors[i])

    ax.set_xlabel(f"PC{dim2+1}")
    ax.set_ylabel("Density")
    plt.savefig("pca_density_y.png", dpi=300)

    # Create density plots for the y-axis
    fig, ax = plt.subplots(figsize=(8, 6))

    for i, name in enumerate(nameLabels):
        mask = labels == name  # Boolean mask
        selected_rows = PCAA[u_norm][mask]  # Select rows where labels match

        if selected_rows.shape[0] > 0:  # Ensure there are matching rows
            group_points = selected_rows[:, [dim1, dim2]]  # Select dim1, dim2
            ax.scatter(
                group_points[:, 0],
                group_points[:, 1],
                label=name,
                color=colors[i],
                edgecolors="black",
            )
            if dis == "y":
                sample_names_array = np.array(
                    sample_names
                )  # Ensure sample_names is a NumPy array
                filtered_names = sample_names_array[mask].tolist()
                for j in range(group_points.shape[0]):
                    ax.text(
                        group_points[j, 0],
                        group_points[j, 1],
                        filtered_names[j],
                        fontsize=9,
                        ha="right",
                        va="bottom",
                        color="black",
                    )

        else:
            print(f"Warning: No data points found for group {name}")

    ax.set_xlabel(f"PC{dim1+1}")
    ax.set_ylabel(f"PC{dim2+1}")
    ax.legend()
    plt.savefig("pca_scatter.png")  # Save the output as PNG
    plt.show()

    num_pcs = mw_PCA.shape[2]
    pcs = np.arange(1, num_pcs + 1)  # Principal components (1-based index)

    evaluator_values = mw_PCA[u_norm, 0, :]

    # Identify the min/max indices based on ranking method
    if u_rank == "p":
        best_index = np.argmin(evaluator_values)  # Minimum for p-value
        highlight_color = "skyblue"
        threshold_line = 0.05
    elif u_rank in ["auc", "aupr"]:
        best_index = np.argmax(evaluator_values)  # Maximum for AUC/AUPR
        highlight_color = "skyblue"
        threshold_line = 0.5
    elif u_rank in ["pc", "sc"]:
        abs_values = np.abs(evaluator_values)
        best_index = np.argmax(abs_values)  # Max absolute for correlations
        highlight_color = "black"

    bar_colors = ["steelblue"] * num_pcs
    bar_colors[best_index] = highlight_color

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 12))

    # First subplot: Evaluator bar plot
    sns.barplot(x=pcs, y=evaluator_values, palette=bar_colors, ax=axes[0])
    axes[0].set_xlabel("Principal Component")
    axes[0].set_ylabel(evaluat)
    axes[0].set_title(
        f"{evaluat} over Principal Components in {ttl} for norm: {u_norm_n}"
    )
    axes[0].axhline(threshold_line, color="red", linestyle="dashed", linewidth=1)

    # Annotate the best-performing PC
    axes[0].annotate(
        f"{evaluator_values[best_index]:.2f}",
        xy=(best_index + 1, evaluator_values[best_index] + 0.02),
        color="red",
        fontsize=10,
        ha="center",
    )

    # Second subplot: Explained variance bar plot
    explained_variance = explained[u_norm]
    sns.barplot(x=pcs, y=explained_variance, color="steelblue", ax=axes[1])
    axes[1].set_xlabel("Principal Component")
    axes[1].set_ylabel("Explained Variance (%)")
    axes[1].set_title("Explained Variance for the Respective Principal Components")

    plt.tight_layout()
    plt.savefig("bar_plots.png", dpi=300)
    plt.show()

    cut_off = get_user_cutoff()

    a = np.zeros(numbLabels)

    # Compute median PCA values for each label
    for i, name in enumerate(nameLabels):
        a[i] = np.median(PCAA[u_norm][labels == name, dim1])

    # Get indices that would sort 'a'
    I = np.argsort(a)

    # Initialize color list
    col = [None] * numbLabels
    col[I[0]] = default_colors[0]  # Smallest median
    col[I[-1]] = default_colors[1]  # Largest median

    Edges = []
    Nodes = []
    pc_corr = []
    x2 = []
    cutoff_f = []

    for i, c in enumerate(cut_off):
        pc[u_norm] = pc[u_norm].T
        print(f"{(u_dim - 1)=}")
        pc_corr_res = C_corr(
            norm_matrices[u_norm], pc[u_norm][:, u_dim - 1], feat_names, c
        )  # Python is zero-indexing so need o get the correct PC loadings column

        Edges.append({"cutoff": c, "data": pc_corr_res["Edges"]})
        pc_corr.append({"cutoff": c, "data": pc_corr_res["pc_corr"]})
        x2.append({"cutoff": c, "data": pc_corr_res["x1"]})
        cutoff_f.append(pc_corr_res["cutoff_f"])

        samplCol = match_V_samplCol(
            col, pc_corr_res["x1"], labels, nameLabels, pc_corr_res["Nodes"]
        )

        Nodes.append(
            {
                "cutoff": c,
                "data": pd.DataFrame(
                    {
                        "Node": pc_corr_res["Nodes"]["Node"],
                        "Colour": samplCol["NodeColor"],
                        "Loading (V)": pc_corr_res["Nodes"]["Loading (V)"],
                    }
                ),
            }
        )

        print(
            f"\nAt cut-off {round(c, 2)}, the PC-corr network has {Nodes[i]['data'].shape[0]} nodes and {Edges[i]['data'].shape[0]} edges.\n"
        )

    save_PC_corr_results_to_excel(cut_off, Edges, Nodes)

    hide_negative_links = get_user_hide_negative_links()

    visualize_network(Edges, Nodes, hide_negative_links)

    # trustworthiness = get_user_trustworthiness(u_lab)

    # if trustworthiness:

    # if u_lab == "c":
    #     numb_rand = 1000  # Number of random permutations
    #     labels_rp = [np.random.permutation(labels) for _ in range(numb_rand)]

    #     trust_cols = ["Trustworthiness_p", "Trustworthiness_auc", "Trustworthiness_aupr"]
    #     for col in trust_cols:
    #         if col not in results.columns:
    #             results[col] = ""

    #     for u_norm_n in results.iloc[:, 5].unique():  # Column 5 contains normalization methods
    #         for u_cent_n in results.iloc[:, 6].unique():  # Column 6 contains centering options
    #             for u_dim in results.iloc[:, 7].unique():  # Column 7 contains dimensions

    #                 print(f"{u_norm_n}")
    #                 print(f"{u_cent_n}")
    #                 print(f"{u_dim}")
    #                 print(f"{u_norm}")

    #                 chosen_opt = (
    #                     (results.iloc[:, 5] == u_norm_n)
    #                     & (results.iloc[:, 6] == u_cent_n)
    #                     & (results.iloc[:, 7] == u_dim)
    #                 )
    #                 idx_opt = np.where(chosen_opt)[0]

    #                 if len(idx_opt) == 0:
    #                     continue  # Skip if no matching rows

    #                 mw_opt = results.iloc[idx_opt, 0].values
    #                 auc_opt = results.iloc[idx_opt, 1].values
    #                 aupr_opt = results.iloc[idx_opt, 2].values

    #                 # Compute trustworthiness values for this combination
    #                 pVal_opt = significance_segregation(
    #                     numb_rand, PCAA[u_norm], u_dim-1, labels_rp, nameLabels, numbLabels,
    #                     [mw_opt, auc_opt, aupr_opt], "all", u_aupr, u_aupr_r
    #                 )

    #                 pVal_mw_opt, pVal_auc_opt, pVal_aupr_opt = pVal_opt

    #                 # Get the correct column indices
    #                 idx_p, idx_auc, idx_aupr = [results.columns.get_loc(col) for col in trust_cols]

    #                 # Update results with computed p-values
    #                 results.iloc[idx_opt, idx_p] = pVal_mw_opt #.astype(str)
    #                 results.iloc[idx_opt, idx_auc] = pVal_auc_opt #.astype(str)
    #                 results.iloc[idx_opt, idx_aupr] = pVal_aupr_opt #.astype(str)

    #     filename = "results.xlsx"
    #     with pd.ExcelWriter(filename, engine="openpyxl", mode="w") as writer:
    #         results.to_excel(writer, sheet_name="PCA results", index=False)

    # print("Results updated with trustworthiness results and saved to 'results.xlsx'.")
    return Edges, Nodes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run PC-Corr analysis from a CSV file."
    )
    parser.add_argument("csv_file", type=str, help="Path to the input CSV file.")
    parser.add_argument(
        "--sample_labels",
        type=str,
        default=None,
        nargs="+",
        help="List of sample labels (space-separated).",
    )
    parser.add_argument(
        "--default_colors",
        type=str,
        default="black,red",
        help="Comma-separated string of colors for distinguishing sample groups (default: 'black,red').",
    )

    args = parser.parse_args()

    default_colors_list = args.default_colors.split(",")

    try:
        df = pd.read_csv(args.csv_file)
        print(f"Loaded CSV file: {args.csv_file}")
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        sys.exit(1)

    if not args.sample_labels:
        sample_names = df.iloc[:, 0].tolist()
        sample_labels: list[str | int] = [
            "East Asian" if x in ["CHN", "KOR", "JPN"] else "Non-East Asian"
            for x in sample_names
        ]
    else:
        sample_labels = args.sample_labels

    PC_corr_v2(df, sample_labels, default_colors=default_colors_list)
