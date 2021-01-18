

import matplotlib.pyplot as plt
import itertools
import numpy as np
from joblib import Parallel, delayed
from svbench.io_tools import CallSet, quantify
from sys import stderr
import numpy as np

__all__ = ["score", "reference_calls_found", "plot", ]


def score(ref_data, query_data, rescore=True, force_intersection=False, reciprocal_overlap=0., stratify=False, n_jobs=1,
          ref_size_bins=(30, 50, 500, 5000, 260000000), allow_duplicate_tp=True, pct_size=0.2):

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
        query = quantify(targets[query.dataset], query, force_intersection, reciprocal_overlap, stratify=stratify,
                         ref_size_bins=ref_size_bins, allow_duplicate_tp=allow_duplicate_tp, pct_size=pct_size)
        print("-" * 45, file=stderr)
        return query

    elif isinstance(query_data, list) or isinstance(query_data, dict):

        if isinstance(query_data, list):
            vals = query_data
        else:
            vals = query_data.values()

        jobs = []
        for query in vals:
            print(f"Scoring dataset={query.dataset} caller={query.caller}", file=stderr)
            if query.dataset not in targets:
                raise ValueError("Reference data not found for ", query.dataset)

            if not rescore:
                continue
            if query.scores is not None:
                # Todo Dont need to use interval tree, just use labelled df
                pass

            if n_jobs == 1:
                print(f"Score table caller={query.caller} against dataset={query.dataset}", file=stderr)
                quantify(targets[query.dataset], query, force_intersection, reciprocal_overlap, stratify=stratify,
                         ref_size_bins=ref_size_bins, allow_duplicate_tp=allow_duplicate_tp, pct_size=pct_size)
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
    if len(table) == 0:
        return None
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
    elif v == "Duplication":
        sdtp = np.in1d(table["TP"], True).sum()
        if sdtp > 0:
            return np.in1d(table["DTP"], True).sum() / np.in1d(table["TP"], True).sum()
        return 0
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


def plot(query_data, x="TP", y="Precision", y2=None, xlim=None, ylim=None, y2lim=None, show=True, refs=None, save_name=None,
         duplicate_tp=False):
    choices = {'Total', 'TP', 'FP', 'FN', 'Duplication', 'Precision', 'Sensitivity', 'DTP', "F1", "Recall"}
    if x not in choices or y not in choices:
        raise ValueError("x and y must be one of: ", choices)

    if isinstance(query_data, dict):
        query_data = list(query_data.values())

    if any(q.false_negative_indexes is None for q in query_data):
        raise ValueError("quant_tools.score function must be called prior to plotting")

    markers = itertools.cycle(['o', '>', '<', '+', 'D', 's', 11, 10, '*', 'X', 3, '_'])

    plots = {}
    for ref_name, grp in itertools.groupby(sorted(query_data, key=lambda qq: qq.dataset), key=lambda q: q.dataset):

        if refs is not None and ref_name not in refs:
            continue

        grp = list(grp)

        plt.figure(figsize=(7, 4))
        ax = plt.subplot(111)
        plt.subplots_adjust(right=0.6)
        plt.title(ref_name)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax2 = None
        if y2 is not None:
            ax2 = ax.twinx()
        for cs in grp:
            if cs.breaks_df is not None:
                marker = next(markers)
                bed = cs.breaks_df[cs.breaks_df["quantified"]]
                # if not duplicate_tp and "DTP" in bed.columns:
                #     bed = bed[~bed["DTP"]]
                if "strata" not in bed.columns or cs.stratify_range is None:
                    ax.scatter(cs.scores[x], cs.scores[y], label=cs.caller, marker=marker)
                    if ax2:
                        ax.scatter(cs.scores[x], cs.scores[y])
                else:
                    x_val = []
                    y_val = []
                    y_val2 = []
                    for threshold in cs.stratify_range:
                        # df = bed[(bed["strata"] >= threshold) & (~bed["DTP"])]
                        df = bed[bed["strata"] >= threshold]
                        x_val.append(calc_score(df, x))
                        y_val.append(calc_score(df, y))
                        if ax2:
                            y_val2.append(calc_score(df, y2))

                    # ax.scatter(x_val, y_val, alpha=0.6, marker=marker, s=20)
                    # ax.plot(x_val, y_val, label=cs.caller, alpha=0.6)
                    ax.plot(x_val, y_val, label=cs.caller, marker=marker, markersize=4, alpha=0.6)

                    if ax2:
                        # ax2.scatter(x_val, y_val2, alpha=0.6, marker=marker, s=20)
                        ax2.plot(x_val, y_val2, linestyle='dashed', label=cs.caller, marker=marker, markersize=4, alpha=0.6)
            else:
                next(markers)

        ax.set_xlabel(x)
        ax.set_ylabel(y)
        if ax2:
            ax2.set_ylabel(y2)
        if xlim is not None:
            plt.xlim(*xlim)
        if ylim is not None:
            ax.set_ylim(*ylim)
        if y2lim is not None:
            ax2.set_ylim(*y2lim)
        # if ylim is not None:
        #     plt.ylim(*ylim)
        plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
        plots[ref_name] = plt
        if save_name:
            plt.savefig(save_name)
    if show:
        plt.show()
    return plots


if __name__ == "__main__":

    pass
