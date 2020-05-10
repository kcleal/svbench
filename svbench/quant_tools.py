

import matplotlib.pyplot as plt
import itertools
import numpy as np
from joblib import Parallel, delayed
from svbench.io_tools import CallSet, quantify
from sys import stderr

__all__ = ["score", "reference_calls_found", "plot", ]


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
