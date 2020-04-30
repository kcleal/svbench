
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from svbench.io_tools import CallSet
from sys import stderr

__all__ = ["score", "reference_calls_found", "plot", ]


def sv_key(chrom, start, chrom2, end):
    # Return a sorted record
    if chrom2 < chrom:
            return chrom2, end, chrom, start
    if end < start:
            return chrom2, end, chrom, start
    else:
        return chrom, start, chrom2, end


def intersecter(tree, chrom, start, end):
    if chrom in tree:
        overlaps = list(tree[chrom].ncls.find_overlap(start, end))
        if len(overlaps) > 0:
            return overlaps


def best_index(out_edges, key="q"):
    max_w = 0
    min_dis = 0
    best_i = 0
    for index, (u, v, d) in enumerate(out_edges):
        if d["weight"] is None or d["weight"] > max_w:
            max_w = d["weight"]
            min_dis = d["dis"]
            best_i = index

        elif d["weight"] == max_w and d["dis"] < min_dis:
            max_w = d["weight"]
            min_dis = d["dis"]
            best_i = index

    m = out_edges[best_i]
    if m[0][0] == key:  # Get query index
        irow = {m[0][1]: m[1][1]}
    else:
        irow = {m[1][1]: m[0][1]}

    others = {}
    for i, (u, v, d) in enumerate(out_edges):
        if i != best_i:
            if u[0] == key:
                others[u[1]] = v[1]
            else:
                others[v[1]] = u[1]

    return irow, others


def quantify(ref_data, data, force_intersection=True, reciprocal_overlap=0., show_table=True, stratify=False):

    # Build a Maximum Bipartite Matching graph:
    # https://www.geeksforgeeks.org/maximum-bipartite-matching/

    G = nx.Graph()

    ref_bedpe = ref_data.breaks_df
    dta = data.breaks_df

    tree = ref_data.tree
    if tree is None:
        raise ValueError("Interval tree has not been created, call add_intervals first")

    # # Link query calls to reference calls
    for query_idx, chrom, start, chrom2, end, w in zip(dta.index, dta["chrom"], dta["start"], dta["chrom2"], dta["end"],
                                                       dta["w"]):

        chrom, start, chrom2, end = sv_key(chrom, start, chrom2, end)

        ol_start = intersecter(tree, chrom, start, start + 1)
        if not ol_start:
            continue

        ol_end = intersecter(tree, chrom2, end, end + 1)
        if not ol_end:
            continue

        # Get the ref_data index
        common_idxs = set([i[2] for i in ol_start]).intersection([i[2] for i in ol_end])
        if len(common_idxs) == 0:
            continue

        # Choose an index by highest weight/lowest total distance, meeting reciprocal_overlap threshold
        min_d = 1e12
        chosen_index = None
        for index in common_idxs:
            ref_row = ref_bedpe.loc[index]
            ref_chrom, ref_start, ref_chrom2, ref_end = sv_key(ref_row["chrom"], ref_row["start"], ref_row["chrom2"],
                                                               ref_row["end"])

            # Make sure chromosomes match
            if chrom != ref_chrom or chrom2 != ref_chrom2:
                continue

            # If intra-chromosomal, check reciprocal overlap
            if chrom == chrom2:
                query_size = end - start + 1e-3
                ref_size = ref_end - ref_start + 1e-3

                ol = float(max(0, min(end, ref_end) - max(start, ref_start)))

                if (ol / query_size < reciprocal_overlap) or (ol / ref_size < reciprocal_overlap):
                    continue

                if force_intersection and ol == 0:
                    continue

            dis = abs(ref_start - start) + abs(ref_end - end)
            if dis < min_d:
                min_d = dis
                chosen_index = index

        if chosen_index is not None:
            G.add_edge(('t', chosen_index), ('q', query_idx), dis=min_d, weight=w)

    good_idxs = {}
    duplicate_idxs = {}
    # Make partitions bipartite-matching i.e. one reference call matched to one query call
    for sub in nx.connected_components(G):
        sub = list(sub)

        if len(sub) == 2:  # One ref matches one query, easy case
            if sub[0][0] == "q":
                good_idxs[sub[0][1]] = sub[1][1]
            else:
                good_idxs[sub[1][1]] = sub[0][1]
            continue

        else:
            bi_count = Counter([i[0] for i in sub])
            if bi_count["t"] == 1 or bi_count["q"] == 1:  # Choose best edge
                ref_node = [i for i in sub if i[0] == "t"][0]
                out_edges = list(G.edges(ref_node, data=True))
                found, others = best_index(out_edges, key="q")
                good_idxs.update(found)
                duplicate_idxs.update(others)

            else:  # Do multi with the maximum-bipartite-matching
                # need to do "q" == 1 and "t" == 2, then  algorithm DELLY
                print("Maximum bipartite matching not implemented", file=stderr)
                quit()
        continue

    missing_ref_indexes = set(ref_bedpe.index).difference(set(good_idxs.values()))

    duplicate_idxs = {k: v for k, v in duplicate_idxs.items() if k not in good_idxs}
    index = data.breaks_df.index  # Use the original IDs?
    data.breaks_df["TP"] = [i in good_idxs for i in index]
    data.breaks_df["DTP"] = [i in duplicate_idxs for i in index]

    data.breaks_df["ref_index"] = [good_idxs[i] if i in good_idxs else duplicate_idxs[i] if i in duplicate_idxs
                                  else None for i in index]
    data.breaks_df["FP"] = [not i and not j for i, j in zip(data.breaks_df["DTP"], data.breaks_df["TP"])]

    if stratify:
        rng = data.stratify_range
    else:
        rng = [None]

    dta = data.breaks_df
    ts = []
    for threshold in rng:

        if threshold is None:
            t = {"Total": len(dta),
                 "Ref": len(ref_bedpe),
                 "DTP": len(duplicate_idxs),
                 "TP": len(good_idxs),
                 "FP": len(index) - len(good_idxs) - len(duplicate_idxs),
                 "FN": len(missing_ref_indexes),
                 "T >=": None
                 }
            sub_total = t["Total"] - t["DTP"]
            t.update({"Precision": round(float(t["TP"]) / sub_total, 4),  # Note DTP are not included
                      "Sensitivity": round(float(t["TP"]) / t["Ref"], 4),
                      "Recall": round(float(t["TP"]) / t["Ref"], 4)})
            t.update({"F1": round(2 * ((t["Precision"] * t["Recall"]) / (t["Precision"] + t["Recall"])), 4)})
            ts.append(t)

        else:


            df = dta[dta["strata"] >= threshold]

            t = {"Total": len(df),
                 "Ref": len(ref_bedpe),
                 "DTP": np.sum(np.in1d(df["DTP"], True)),
                 "TP": np.sum(np.in1d(df["TP"], True)),
                 "FP": np.sum(np.in1d(df["FP"], True)),
                 "FN": None,
                 "T >=": threshold
                 }
            sub_total = t["Total"] - t["DTP"]
            t.update({"Precision": round(float(t["TP"]) / sub_total, 4),
                      "Recall": round(float(t["TP"]) / t["Ref"], 4)})
            t.update({"F1": round(2 * ((t["Precision"] * t["Recall"]) / (t["Precision"] + t["Recall"])), 4)})
            ts.append(t)

    data.scores = pd.DataFrame.from_records(ts)[["T >=", "Ref", "Total", "TP", "FP", "DTP", "FN", "Precision",
                                                 "Recall", "F1"]]
    if show_table:
        try:
            print(data.scores.to_markdown(), file=stderr)
        except:
            print(data.scores.to_string(), file=stderr)
    data.false_negative_indexes = missing_ref_indexes

    return data


def score(ref_data, query_data, rescore=True, force_intersection=True, reciprocal_overlap=0., stratify=False, n_jobs=1):

    if isinstance(ref_data, CallSet):
        targets = {ref_data.dataset: ref_data}
    elif isinstance(ref_data, list):
        targets = {v.dataset: v for v in ref_data}
    elif isinstance(ref_data, dict):
        targets = ref_data
    else:
        raise ValueError("ref_data argument must be of type CallSet, list or dict")

    if isinstance(query_data, CallSet):
        query = query_data
        print(f"Score table caller={query.caller} against dataset={query.dataset}", file=stderr)
        query = quantify(targets[query.dataset], query, force_intersection, reciprocal_overlap, stratify=stratify)
        print("-" * 45, file=stderr)
        return query

    elif isinstance(query_data, list) or isinstance(query_data, dict):

        if isinstance(query_data, list):
            vals = query_data
        else:
            vals = query_data.values()

        jobs = []
        for query in vals:

            if query.dataset not in targets:
                raise ValueError("Reference data not found for ", query.dataset)

            if not rescore:
                continue
            if query.scores is not None:
                # Todo Dont need to use interval tree, just use labelled df
                pass

            if n_jobs == 1:
                print(f"Score table caller={query.caller} against dataset={query.dataset}", file=stderr)
                quantify(targets[query.dataset], query, force_intersection, reciprocal_overlap, stratify=stratify)
                print("-"*45, file=stderr)

            else:
                jobs.append((targets[query.dataset], query, force_intersection, reciprocal_overlap, False, stratify))

        if jobs:
            vals = Parallel(n_jobs=n_jobs)(delayed(quantify)(*args) for args in jobs)

        if isinstance(query_data, list):
            return vals
        else:
            return {k: v for k, v in zip(query_data.keys(), vals)}

    else:
        raise ValueError("query_data argument must be of type CallSet, list or dict")


def calc_score(table, v, ref=None):
    # Assume pandas DataFrame
    if v == "Total":
        return len(table)
    elif v in ("TP", "FP", "DTP"):
        return np.in1d(table[v], True).sum()
    elif v == "Precision":
        return np.in1d(table["TP"], True).sum() / float(len(table))
    elif v == "Sensitivity":
        if ref is None:
            raise ValueError("ref data must be provided to calculate sensitivity")
        return table["TP"] / len(ref)
    else:
        raise ValueError("Not implemented for table: {}".format(v))


def reference_calls_found(ref_data, query_data):
    """
    Add columns to reference DataFrame showing which calls were identified/missed in query call sets
    :param ref_data: Reference SV calls, CallSet
    :param query_data: Query_data calls, CallSet
    :return: Ref_data, CallSet
    """
    if not isinstance(ref_data, CallSet) or not isinstance(query_data, CallSet):
        raise ValueError("ref_data and query_data must be an instance of CallSet, not {}, {}".format(type(ref_data),
                                                                                                     type(query_data)))
    if ref_data.dataset != query_data.dataset:
        raise ValueError("dataset attribute must match for ref_data and query_data")

    found_idxs = set(query_data.breaks_df["ref_index"].dropna())
    ref_data.breaks_df[query_data.caller] = [1 if i in found_idxs else 0 for i in ref_data.breaks_df.index]
    return ref_data


def plot(query_data, x="TP", y="Precision", xlim=None, ylim=None, show=True, refs=None, save_prefix=None,
         duplicate_tp=False):
    choices = {'Total', 'TP', 'FP', 'FN', 'Precision', 'Sensitivity', 'DTP', "F1", "Recall"}
    if x not in choices or y not in choices:
        raise ValueError("x and y must be one of: ", choices)

    if isinstance(query_data, dict):
        query_data = list(query_data.values())

    if any(q.false_negative_indexes is None for q in query_data):
        raise ValueError("quant_tools.score function must be called prior to plotting")

    markers = itertools.cycle(['o', '>', '<', '+', 'D', 's', 11, 10, '*', 'X', 3, '_'])

    plots = {}
    for ref_name, grp in itertools.groupby(sorted(query_data, key=lambda qq: (qq.dataset, qq.caller)),
                                           key=lambda q: q.dataset):

        if refs is not None and ref_name not in refs:
            continue

        grp = list(grp)

        plt.figure(figsize=(7, 4))
        ax = plt.subplot(111)
        plt.subplots_adjust(right=0.6)
        plt.title(ref_name)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        for cs in grp:
            bed = cs.breaks_df
            if not duplicate_tp and "DTP" in bed.columns:
                bed = bed[~bed["DTP"]]
            if "strata" not in bed.columns:
                ax.scatter(cs.scores[x], cs.scores[y], label=cs.caller, marker=next(markers))

            else:
                x_val = []
                y_val = []
                for threshold in cs.stratify_range:
                    # df = bed[(bed["strata"] >= threshold) & (~bed["DTP"])]
                    df = bed[bed["strata"] >= threshold]
                    x_val.append(calc_score(df, x))
                    y_val.append(calc_score(df, y))

                ax.scatter(x_val, y_val, alpha=0.6, marker=next(markers), s=20)
                ax.plot(x_val, y_val, label=cs.caller, alpha=0.6)

        plt.xlabel(x)
        plt.ylabel(y)
        if xlim is not None:
            plt.xlim(*xlim)
        if ylim is not None:
            plt.ylim(*ylim)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plots[ref_name] = plt
        if save_prefix:
            plt.savefig(save_prefix + str(ref_name) + ".pdf")
    if show:
        plt.show()
    return plots


if __name__ == "__main__":

    pass
