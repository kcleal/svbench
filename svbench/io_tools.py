import click
import pandas as pd
import ncls
from collections import defaultdict
import numpy as np
import vcf
import operator
import datetime
import pkg_resources
import copy
import pickle
from sys import stderr, stdout, stdin
import gzip
import resource
import time
from io import StringIO


__all__ = ["Col", "CallSet", "concat_dfs"]


class Col:
    """This is a column parser class. The input column must be index-able by 'col' argument. Subfields are accessed
    using the 'key' argument. Values can be further encoded by providing a dict with the required mappings using the
    'encoding' argument. The 'bins' argument can be used to stratify input into predefined bins.
    :param col: The key of the primary data field
    :rtype col: str, optional
    :param key: The key of the secondary data field, optional
    :rtype key: str, optional

    """
    def __init__(self, col, key=None, encoding=None, bins=None, norm=None, op=None, thresh=None):
        self.col = col
        self.key = key
        self.encoding = encoding
        self.bins = bins
        self.norm = norm
        self.op = op
        self.thresh = thresh

    def __repr__(self):
        return f"svbench.Col(col={self.col}, key={self.key}, encoding={self.encoding}, bins={self.bins}, norm={self.norm}, op={self.op}, thresh={self.thresh})"


class Operate:
    def __init__(self):
        self.opps = {k: eval('operator.' + k) for k in dir(operator) if "_" not in k}

    def test(self, o, a, b):
        if isinstance(a, list):
            if len(a) == 0:
                a = None
            else:
                a = a[0]
        if isinstance(b, list):
            if len(b) == 0:
                a = None
            else:
                b = b[0]
        if isinstance(a, str) and a.isnumeric():
            a = float(a)
        if isinstance(b, str) and b.isnumeric():
            b = float(b)
        return self.opps[o](a, b)


class NSV:
    def __init__(self, starts, ends, ids):
        self.starts = starts
        self.ends = ends
        self.ids = ids
        self.ncls = ncls.NCLS(starts, ends, ids)

    def __reduce__(self):  # Problems pickling original ncls object, so made a new reduce method
        return self.__class__, (self.starts, self.ends, self.ids)


def get_interval_arrays(regions, slop, bedpe=False):
    chrom_interval_start = defaultdict(list)
    chrom_interval_end = defaultdict(list)
    chrom_interval_index = defaultdict(list)

    # Make a region for each breakpoint
    if not bedpe:
        for c1, s, c2, e, index in regions:
            chrom_interval_start[c1].append(s - slop)
            chrom_interval_end[c1].append(s + slop)
            chrom_interval_index[c1].append(index)

            chrom_interval_start[c2].append(e - slop)
            chrom_interval_end[c2].append(e + slop)
            chrom_interval_index[c2].append(index)

    else:
        for c1, s1, e1, c2, s2, e2, index in regions:
            chrom_interval_start[c1].append(s1)
            chrom_interval_end[c1].append(e1)
            chrom_interval_index[c1].append(index)

            chrom_interval_start[c2].append(s2)
            chrom_interval_end[c2].append(e2)
            chrom_interval_index[c2].append(index)

    return ({k: np.array(v).astype(int) for k, v in chrom_interval_start.items()},
           {k: np.array(v).astype(int) for k, v in chrom_interval_end.items()},
           {k: np.array(v).astype(int) for k, v in chrom_interval_index.items()})


def make_ncls_table(chrom_interval_start, chrom_interval_end, chrom_interval_index):
    return {k: NSV(chrom_interval_start[k], chrom_interval_end[k], chrom_interval_index[k])
            for k in chrom_interval_start}


def make_tree(regions, slop=250, bedpe=False):
    chrom_interval_start, chrom_interval_end, chrom_interval_index = get_interval_arrays(regions, slop, bedpe)
    return make_ncls_table(chrom_interval_start, chrom_interval_end, chrom_interval_index)


def col_parser(r, col, key, first=True):
    if col == "FORMAT":
        col = r.__getattribute__("samples")
        if len(col) > 1:
            raise ValueError("SVBench only supports one FORMAT column per vcf file, {} found", len(col))
        try:
            ck = col[0][key]
        except IndexError:
            raise ValueError("The FORMAT column is missing a name i.e. #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT   {name}???")

    else:
        try:
            col = r.__getattribute__(col)
        except KeyError:
            raise ValueError("Column argument not understood, did you mean 'INFO'?")

        if key is not None:
            if key in col:
                ck = col[key]
            else:
                ck = None

        else:
            ck = col
        if isinstance(ck, list):
            if len(ck) == 1:  # Get first item of list or use None
                return ck[0]
            elif len(ck) > 1:
                if first:
                    ck = ck[0]
                else:
                    return ck
            else:
                ck = None
    return ck


def check_passed(operations, r, keep):
    passed = True
    for item in keep:
        if item.op is None:
            raise ValueError("Col.op must be set using 'keep' argument e.g. "
                             "Col('INFO', 'SU', op=eq, thresh=4)")
        ck = col_parser(r, item.col, item.key)
        if not operations.test(item.op, ck, item.thresh):
            passed = False
            break
    return passed


def get_strata(r, strat):

    cp = col_parser(r, strat.col, strat.key)
    if cp is None:
        cp = 0
    return float(cp)


def parse_cols(r, item):

    encoding = item.encoding
    col = item.col
    key = item.key

    if encoding:

        if key is not None:
            p = col_parser(r, item.col, item.key)
            if p in encoding:
                return f"{col}:{key}", encoding[p]
            elif None in encoding:
                return f"{col}:{key}", encoding[None]
            else:
                raise ValueError("{}:{} could not be found in encoding".format(col, key), item)

        else:
            if type(r) == vcf.model._Record:
                p = r.__getattribute__(col)
                if isinstance(p, list):
                    if len(p) == 0:
                        p = None
                    else:
                        p = p[0]
            else:
                p = r[col]

            if p in encoding:
                return col, encoding[p]
            else:
                return col, encoding[None]

    else:
        if key is not None:
            p = col_parser(r, col, key, first=False)
            if isinstance(p, list):
                d = dict()
                for index, v in enumerate(p):
                    d[f"{col}:{key}:{index}"] = v
                return d
            else:
                return f"{col}:{key}", p
        else:
            return col, r[col]


def parse_cols_list(r, parse_list, data):

    for item in parse_list:
        d = parse_cols(r, item)
        if isinstance(d, dict):
            data.update(d)
        else:
            data[d[0]] = d[1]

    return data


def check_args(stratify, weight_field, keep, other_cols, size_range):

    if size_range is not None:
        if size_range[0] is not None and size_range[1] is not None:
            assert size_range[0] < size_range[1]

    for item in (keep, weight_field, stratify, other_cols):
        if item is None:
            continue
        if isinstance(item, list):
            for item2 in item:
                assert isinstance(item2, Col)
        else:
            assert isinstance(item, Col)

    if stratify is not None and (stratify.bins is None or not hasattr(stratify.bins, '__iter__')):
        raise ValueError(
            "Stratify must have an iterable for 'bins' argument e.g. Col('FORMAT', 'DV', bins=range(0, 60, 5))")


def filter_by_size(df, size_range=(None, None)):
    l_before = len(df)
    if size_range[0] is not None or size_range[1] is not None:
        min_size, max_size = size_range
        drop = []
        for idx, start, end, chrom1, chrom2 in zip(df.index, df["start"], df["end"], df["chrom"], df["chrom2"]):
            if chrom1 == chrom2:
                s = abs(end - start)
                if min_size is not None and s < min_size:
                    drop.append(idx)
                elif max_size is not None and s > max_size:
                    drop.append(idx)
        print("Length before/after filter: {}, {}".format(l_before, l_before - len(drop)), file=stderr)
        return df.drop(index=drop, inplace=False)
    return df


class CallSet:
    def __init__(self, dataset=None, caller=None, **kwargs):

        # Primary attributes
        self.dataset = dataset
        self.caller = caller
        self.name = None
        self.kwargs = kwargs
        self.required = {}  # Set any number of required arguments. Note, arguments are accessed via kwargs

        # Secondary attributes
        self.bedpe = False
        self.path = None
        self.tree = None
        self.breaks_df = None
        self.extra_cols = None
        self.style = "GIAB"
        self.add_weight = False,
        self.weight_field = None,
        self.no_translocations = True
        self.allowed_svtypes = None
        self.keep = None
        self.stratify = None
        self.stratify_range = None
        self.allowed_chroms = None
        self.min_size = None
        self.max_size = None
        self.other_cols = None
        self.model = None
        self.slop = None
        self.break_cols = "chrA,posA,chrB,posB"
        self.sep = ","
        self.drop_first_row = False
        self.svtype_field = "svtype"
        self.scores = None
        self.false_negative_indexes = None
        self.kind = None
        self.new_col = None
        self.add_to = None

        self.temp = None

        # Track defaults for persistence
        self.default_params = {k: v for k, v in self.__dict__.items() if k in {"bedpe",
                                                                               "style",
                                                                               "add_weight",
                                                                               "weight_field",
                                                                               "no_translocations",
                                                                               "allowed_svtypes",
                                                                               "keep",
                                                                               "stratify",
                                                                               "stratify_range",
                                                                               "allowed_chroms",
                                                                               "min_size",
                                                                               "max_size",
                                                                               "other_cols",
                                                                               "slop",
                                                                               "break_cols",
                                                                               "sep",
                                                                               "drop_first_row",
                                                                               "svtype_field"}}

        self.meta_data = {"scikit-learn": pkg_resources.get_distribution("scikit-learn").version,
                          "svbench": pkg_resources.get_distribution("svbench").version,
                          "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

        # Update any provided kwargs
        if kwargs:
            self.set_properties(kwargs)

    def __add__(self, other):
        if not isinstance(other, CallSet):
            raise ValueError("Calling add on {}, should be instance of CallSet".format(type(other)))

        self.breaks_df["source"] = [self.caller] * len(self.breaks_df)
        other_df = other.breaks_df.copy(deep=True)
        other_df["source"] = [other.caller] * len(other_df)
        self.breaks_df = pd.concat([self.breaks_df, other_df])
        self.tree = None
        return self

    def __sub__(self, other):
        if not isinstance(other, CallSet):
            raise ValueError("Calling add on {}, should be instance of CallSet".format(type(other)))
        if self.tree is None or other.tree is None:
            raise ValueError("Call add_intervals on both CallSet's before subtracting")

        this_df = self.breaks_df
        bad_indexes = set([])
        starts, ends, indexes = get_interval_arrays(list(zip(this_df["chrom"], this_df["start"], this_df["chrom2"],
                                                  this_df["end"], this_df.index)), self.slop)

        other_tree = other.tree
        for chrom in starts.keys():
            if chrom in other_tree:
                l_idxs, r_idxs = other_tree[chrom].ncls.all_overlaps_both(starts[chrom], ends[chrom], indexes[chrom])
                bad_indexes |= set(l_idxs)

        self.breaks_df = this_df.drop(list(bad_indexes))
        self.add_intervals(self.slop)
        return self

    def __iand__(self, other):
        if not isinstance(other, CallSet):
            raise ValueError("Calling add on {}, should be instance of CallSet".format(type(other)))
        if self.tree is None or other.tree is None:
            raise ValueError("Call add_intervals on both CallSet's before subtracting")
        this_df = self.breaks_df
        good_indexes = set([])
        starts, ends, indexes = get_interval_arrays(list(zip(this_df["chrom"], this_df["start"], this_df["chrom2"],
                                                             this_df["end"], this_df.index)), self.slop)

        other_tree = other.tree
        for chrom in starts.keys():
            if chrom in other_tree:
                l_idxs, r_idxs = other_tree[chrom].ncls.all_overlaps_both(starts[chrom], ends[chrom], indexes[chrom])
                good_indexes |= set(l_idxs)

        self.breaks_df = this_df.loc[sorted(good_indexes)]
        self.add_intervals(self.slop)
        return self

    def add_intervals(self, slop=250):
        self.slop = slop
        df = self.breaks_df

        self.tree = make_tree(list(zip(df["chrom"], df["start"], df["chrom2"], df["end"], df.index)), slop=slop)
        return self

    def reset_arguments(self):
        for k, v in self.default_params.items():
            self.__setattr__(k, v)

    def set_args(self, args):
        # Update existing class attributes with any new ones
        for k, v in args.items():
            if k != "self":
                a = self.__getattribute__(k)
                if a is None and v is not None:
                    self.__setattr__(k, v)
                elif v is not None:
                    if k in self.default_params and v != self.default_params[k]:
                        self.__setattr__(k, v)

    def properties(self):
        return self.kwargs

    def required(self):
        return self.required

    def set_properties(self, d):
        unset = set(self.required)
        new_args = set([])
        for k, v in d.items():
            if k in self.__dict__:
                self.__setattr__(k, v)
            else:
                if k in unset:
                    unset.remove(k)
                else:
                    new_args.add(k)
                self.kwargs[k] = v
        if len(unset) > 0:
            raise KeyError(f"All preset key word arguments must be set when calling new. Keys missing: {list(unset)}")
        if len(new_args) > 0:
            print(f"Warning: unexpected arguments {list(new_args)}", file=stderr)
        return self

    def deepcopy(self):
        return copy.deepcopy(self)

    def new(self, **kwargs):
        n = self.deepcopy()
        n.set_properties(kwargs)
        return n

    def filter_by_size(self, min_size=None, max_size=None, inplace=False):
        if inplace:
            cs = self
        else:
            cs = self.deepcopy()
        done = True
        if cs.min_size != min_size:
            cs.min_size = min_size
            done = False
        if cs.max_size != max_size:
            cs.max_size = max_size
            done = False
        if done or (min_size is None and max_size is None):  # Nothing to be done
            return cs

        l_before = len(cs.breaks_df)
        df = cs.breaks_df
        drop = []
        for idx, start, end, chrom1, chrom2 in zip(df.index, df["start"], df["end"], df["chrom"], df["chrom2"]):
            if chrom1 == chrom2:
                s = abs(end - start)
                if min_size is not None and s < min_size:
                    drop.append(idx)
                elif max_size is not None and s > max_size:
                    drop.append(idx)
        print("Filtered by size, caller={}, dataset={} rows before {}, after {}".format(self.caller, self.dataset,
                                                                                        l_before, l_before - len(drop)),
              file=stderr)
        df.drop(index=drop, inplace=True)
        return cs

    def score_table(self):
        if self.scores is not None:
            print(f"Score table caller={self.caller} against dataset={self.dataset}", file=stderr)
            try:
                print(pd.DataFrame.from_records([self.scores], index=None).to_markdown(), file=stderr)
            except:
                print(pd.DataFrame.from_records([self.scores], index=None).to_string(), file=stderr)

    def set_strata(self, stratify_col, stratify_range):

        if not hasattr(stratify_range, '__iter__'):
            raise ValueError("stratify_range must be an iterable")

        self.breaks_df["strata"] = self.breaks_df[stratify_col]
        self.stratify_range = stratify_range
        return self

    def predict_proba(self, model=None, col_name="PROB", show_stats=True, density=False, bins=20, show_model=True):

        if model is None:
            model = self.model
        if model is None:
            raise ValueError("No model to use")
        if show_model:
            print("Model:", model, file=stderr)

        # Sometimes NaN creep in or infinity
        X = np.nan_to_num(self.breaks_df[self.extra_cols].astype(float))
        y_predict = model.predict_proba(X)[:, 1]

        if show_stats:
            hist, bin_edges = np.histogram(y_predict, bins=bins, range=(0, 1), density=density)
            print(
                f"Distribution of probability values, caller={self.caller} dataset={self.dataset} (n={len(y_predict)})"
                , file=stderr)

            pp = pd.DataFrame.from_records([{round(j, 3): i for i, j in zip(hist, bin_edges)}]).set_index(
                pd.Index(["Calls"]))
            try:
                print(pp.to_markdown(), file=stderr)
            except:
                print(pp.to_string(), file=stderr)

        self.breaks_df[col_name] = y_predict
        return self

    def predict(self, model=None, col_name="PROB", show_stats=True, density=False, bins=2, show_model=True):

        if model is None:
            model = self.model
        if model is None:
            raise ValueError("No model to use")
        if show_model:
            print("Model:", model, file=stderr)

        # Sometimes NaN creep in or infinity
        X = np.nan_to_num(self.breaks_df[self.extra_cols].astype(float))
        y_predict = model.predict(X)

        if show_stats:
            hist, bin_edges = np.histogram(y_predict, bins=bins, range=(0, 1), density=density)
            print(
                f"Distribution of probability values, caller={self.caller} dataset={self.dataset} (n={len(y_predict)})"
                , file=stderr)
            pp = pd.DataFrame.from_records([{round(j, 3): i for i, j in zip(hist, bin_edges)}]).set_index(
                pd.Index(["Calls"]))
            try:
                print(pp.to_markdown(), file=stderr)
            except:
                print(pp.to_string(), file=stderr)

        self.breaks_df[col_name] = y_predict
        return self

    def load(self, path):
        if self.kind is None:
            raise ValueError("The file kind must be provided with the argument 'kind'")
        elif self.kind == "vcf":
            return self.load_vcf(path)
        elif self.kind == "csv":
            return self.load_csv(path)
        elif self.kind == "bed" or self.kind == "bedpe":
            return self.load_bed(path)
        else:
            raise ValueError("Unknown file kind {}".format(self.kind))

    def load_vcf(self, path, style="GIAB", add_weight=False, weight_field=None,
                 no_translocations=True, allowed_svtypes=None, keep=None, stratify=None, allowed_chroms=None,
                 min_size=None, max_size=None,
                 other_cols=None):

        self.path = path
        self.kind = "vcf"

        check_args(stratify, weight_field, keep, other_cols, (min_size, max_size))

        self.set_args(locals())  # Overwrite default params

        # Load instance parameters, these may have been set previously
        style = self.style
        add_weight = self.add_weight
        weight_field = self.weight_field
        no_translocations = self.no_translocations
        allowed_svtypes = self.allowed_svtypes
        keep = self.keep
        stratify = self.stratify
        allowed_chroms = self.allowed_chroms
        min_size = self.min_size
        max_size = self.max_size
        other_cols = self.other_cols

        operations = Operate()

        res = []
        if allowed_svtypes is not None:
            allowed_svtypes = set(allowed_svtypes.split(","))

        if allowed_chroms:
            if not isinstance(allowed_chroms, set):
                allowed_chroms = set(allowed_chroms)

        new_cols = []
        new_index = []

        if path in "-stdin":
            temp = StringIO()
            for line in stdin:
                temp.write(line)
            temp.seek(0)
            self.temp = temp  # Used when writing output

            reader = vcf.Reader(fsock=temp)
        else:
            reader = vcf.Reader(filename=path)

        for r in reader:

            chrom = r.CHROM
            if isinstance(chrom, int) or (isinstance(chrom, str) and chrom[0] != "c"):
                chrom = "chr" + str(chrom)

            if allowed_chroms is not None and chrom not in allowed_chroms:
                continue

            start = int(r.POS)

            svtype = r.INFO["SVTYPE"]

            if allowed_svtypes is not None and svtype not in allowed_svtypes:
                continue

            if keep is not None and not check_passed(operations, r, keep):
                continue

            if svtype == "BND":
                if "_2" in r.ID:  # Skip second part of BND record
                    continue
                chrom2 = r.ALT[0].chr
                end = r.ALT[0].pos
                if end is None:
                    end = start + 1
                else:
                    end = int(r.ALT[0].pos)
            else:
                chrom2 = chrom
                if style == "GIAB":
                    end = int(r.INFO["END"])
                else:
                    raise IOError

            if no_translocations and chrom != chrom2:
                continue

            if chrom == chrom2:
                size = abs(end - start)
                if min_size is not None and size < min_size:
                    continue
                if max_size is not None and size > max_size:
                    continue

            if allowed_chroms is not None and chrom2 not in allowed_chroms:
                continue

            w = None

            if add_weight:
                w = parse_cols(r, weight_field)

                if isinstance(w, dict):
                    raise ValueError("Multiple values in weight field, choose a field with a single value")
                wcol, w = w

            d = {"end": end, "start": start, "chrom": chrom, "chrom2": chrom2, "w": w}
            if stratify is not None:
                d["strata"] = get_strata(r, stratify)

            if other_cols:
                parsed = parse_cols_list(r, other_cols, {})
                if not new_cols:
                    new_cols = list(parsed.keys())
                    # print("newcols", new_cols)
                d.update(parsed)

            res.append(d)
            new_index.append(r.ID)

        if len(res) == 0:
            print("WARNING: empty vcf", file=stderr)
            return self

        df = pd.DataFrame.from_records(res)

        # Normalize new columns, might be more than one column added per Col
        if other_cols is not None:
            item_key = {(k.col, k.key): k for k in other_cols if k.norm is not None}
            for k, v in item_key.items():
                for ec in new_cols:
                    ecs = ec.split(":")

                    if ecs[0] == k or (ecs[0], ecs[1]) == k:
                        df[ec] = v.norm(df[ec], self.kwargs)

        # Keep record of IDs
        df["new_index"] = new_index

        # Order df
        base_cols = ["chrom", "start", "chrom2", "end", "w", "strata"]
        df = df[[i for i in base_cols if i in df.columns] + new_cols]

        self.breaks_df = df
        self.extra_cols = new_cols

        if stratify is not None:
            self.stratify_range = stratify.bins

        print(f"dataset={self.dataset}, caller={self.caller}, loaded rows: {len(df)}", file=stderr)

        return self

    def load_bed(self, path, drop_first_row=False, bedpe=False):

        self.path = path
        if bedpe is True:
            self.kind = "bedpe"
        else:
            self.kind = "bed"

        self.set_args(locals())

        drop_first_row = self.drop_first_row

        df = pd.read_csv(path, comment="#", sep="\t", header=None)
        if df[0].iloc[0] == "Chr" or drop_first_row:
            df.drop([0], inplace=True)

        if self.kind == "bedpe":
            df.rename({0: "chrom", 1: "start1", 2: "end1", 3: "chrom2", 4: "start2", 5: "end2"}, axis="columns",
                      inplace=True)
            if df["chrom"].iloc[0][0] != "c":
                df["chrom"] = ["chr" + str(i) for i in df["chrom"]]

            for col in ("start1", "end2", "start2", "end2"):
                df[col] = df[col].astype(int)

        else:
            df.rename({0: "chrom", 1: "start", 2: "end"}, axis="columns", inplace=True)
            if df["chrom"].iloc[0][0] != "c":
                df["chrom"] = ["chr" + str(i) for i in df["chrom"]]
            df["chrom2"] = df["chrom"]
            df["start"] = df["start"].astype(int)
            df["end"] = df["end"].astype(int)

        self.breaks_df = df
        print(f"dataset={self.dataset}, caller={self.caller}, loaded rows: {len(df)}", file=stderr)

        return self

    def load_csv(self, path, break_cols="chrA,posA,chrB,posB", sep=",", add_weight=False, weight_field=None,
                 allowed_svtypes=None, keep=None, svtype_field="svtype", no_translocations=True,
                 allowed_chroms=None, stratify=None,
                 min_size=None, max_size=None,
                 other_cols=None, drop_first_row=False):

        self.path = path
        self.kind = "csv"

        check_args(stratify, weight_field, keep, other_cols, (min_size, max_size))

        self.set_args(locals())

        # Load instance parameters, these may have been set previously
        break_cols = self.break_cols
        sep = self.sep
        add_weight = self.add_weight
        weight_field = self.weight_field
        svtype_field = self.svtype_field
        drop_first_row = self.drop_first_row
        no_translocations = self.no_translocations
        allowed_svtypes = self.allowed_svtypes
        keep = self.keep
        stratify = self.stratify
        allowed_chroms = self.allowed_chroms
        min_size = self.min_size
        max_size = self.max_size
        other_cols = self.other_cols

        operations = Operate()

        df_in = pd.read_csv(path, sep=sep, index_col=None, comment="#")
        if drop_first_row:
            df_in.drop([0], inplace=True)

        if allowed_svtypes is not None:
            allowed_svtypes = set(allowed_svtypes.split(","))
            df_in = df_in[df_in[svtype_field].isin(allowed_svtypes)]

        if keep is not None:
            for item in keep:
                df_in = df_in[[operations.test(item.op, i, item.thresh) for i in df_in[item.col]]]

        assert break_cols.count(",") == 3
        df = pd.DataFrame()
        for n, col in zip(["chrom", "start", "chrom2", "end"], break_cols.split(",")):
            df[n] = df_in[col]

        if allowed_chroms:
            if not isinstance(allowed_chroms, set):
                allowed_chroms = set(allowed_chroms)
            df = df[df["chrom"].isin(allowed_chroms) & df["chrom2"].isin(allowed_chroms)]

        if no_translocations:
            df = df[df["chrom"] == df["chrom2"]]

        if min_size is not None or max_size is not None:
            drop = []
            for idx, start, end, chrom1, chrom2 in zip(df.index, df["start"], df["end"], df["chrom"], df["chrom2"]):
                if chrom1 == chrom2:
                    s = abs(end - start)
                    if min_size is not None and s < min_size:
                        drop.append(idx)
                    elif max_size is not None and s > max_size:
                        drop.append(idx)

            df.drop(index=drop, inplace=True)

        if stratify is not None:
            df["strata"] = df_in[stratify.col].loc[df.index]
            self.stratify_range = stratify.bins

        if add_weight and weight_field is not None:
            df["w"] = df_in[weight_field.col]
        else:
            df["w"] = [None] * len(df)

        if other_cols:

            index = df.index
            for new_col in other_cols:
                if new_col.encoding is not None:
                    enc = new_col.encoding
                    df[new_col.col] = [enc[j] if j in enc else enc[None] for j in df_in[new_col.col].loc[index]]
                else:
                    df[new_col.col] = df_in[new_col.col].loc[index]

                if new_col.norm is not None:
                    df[new_col.col] = new_col.norm(df[new_col.col], self.kwargs)
            self.extra_cols = [k.col for k in other_cols]

        if len(df) == 0:
            print("Warning: empty vcf", file=stderr)
            return self

        self.breaks_df = df
        print(f"dataset={self.dataset}, caller={self.caller}, loaded rows: {len(df)}", file=stderr)
        return self

    def write_to_vcf(self, path, new_col="PROB", add_to="FORMAT"):

        if isinstance(path, str):
            if path[-3:] == ".gz":
                f_out = gzip.open(path, "wb")
            else:
                f_out = open(path, "w")
        else:
            f_out = path  # assume file handle

        def get_f_in(path):
            if self.temp is not None:
                self.temp.seek(0)  # assume io.StringIO
                return self.temp

            if path[-3:] == ".gz":
                return gzip.open(path, "rb")
            return open(path, "r")

        if new_col is not None and not (add_to == "FORMAT" or add_to == "INFO"):
            raise ValueError("add_to argument must be FORMAT or INFO")

        col_vals = {idx: v for idx, v in zip(self.breaks_df.index, self.breaks_df[new_col])}

        replace_val = False
        write_h_line = True
        write_h = True
        header_lines = []
        last_add_index = 0
        position = 0
        with f_out, get_f_in(self.path) as f_in:

            for findex, line in enumerate(f_in):
                if line[0] == "#":

                    if f'##{add_to}' in line:
                        last_add_index = findex
                        if f'##{add_to}=<ID={new_col},' in line:  # Replace line
                            replace_val = True
                            write_h_line = False
                            header_lines.append(f'##{add_to}=<ID={new_col},Number=1,Type=Float,Description="Probability of event, added by SVBench">\n')
                            continue

                    header_lines.append(line)
                    continue

                else:
                    if write_h:
                        if write_h_line:  # this line might already have been replaced, so no need to write again
                            header_lines.insert(last_add_index + 1,
                                f'##{add_to}=<ID={new_col},Number=1,Type=Float,Description="Probability of event, added by SVBench">\n')

                        f_out.writelines(header_lines)
                        write_h = False

                    l = line.split("\t")
                    r_id = position

                    if r_id in col_vals:
                        val = str(col_vals[r_id])
                        if val == "1.0":
                            val = "1"
                    else:
                        val = "0"
                    if add_to == "FORMAT":
                        vidx = None
                        if replace_val:
                            l8s = l[8].split(":")
                            vidx = [idx for idx, vv in enumerate(l8s) if vv == new_col][0]
                            if vidx == len(l8s)-1:
                                vidx = None  # Place at end
                        else:
                            l[8] += ":" + new_col
                        # Value is added to each record in the row
                        for i in range(9, len(l)):
                            if i == len(l) - 1 and vidx is None:
                                l[i] = l[i].strip() + ":" + val + "\n"
                            elif vidx is None:
                                l[i] = l[i] + ":" + val
                            else:
                                lis = l[i].split(":")
                                lis[vidx] = val
                                l[i] = ":".join(lis)

                        f_out.writelines("\t".join(l))
                    elif add_to == "INFO":
                        l[7] += f";{new_col}={val}"
                        f_out.writelines("\t".join(l))

                    position += 1

    def save_as(self, path, new_col=None, add_to=None, format=None):
        if new_col is None and self.new_col is not None:
            new_col = self.new_col
        if add_to is None and self.add_to is not None:
            add_to = self.add_to

        if format is None and self.kind == "vcf":
            self.write_to_vcf(path, new_col, add_to)
        if format is None and self.kind == "csv":
            raise ValueError("Not implemented currently")

        else:
            raise ValueError(self.kind, "not in {csv, vcf}")

    def save_intervals(self, path, format="bed"):
        if format not in {"bed", "bedpe"}:
            raise ValueError("format must be bed or bedpe")
        if self.breaks_df is None:
            raise ValueError("breaks_df object has not been constructed")
        if self.slop is None:
            slop = 1
        else:
            slop = self.slop
        with open(path, "w") as out:

            v = []
            for idx, r in self.breaks_df.iterrows():
                if format == "bed":

                    v.append((r['chrom'], r['start'] - slop, r['start'] + slop, idx))
                    v.append((r['chrom2'], r['end'] - slop, r['end'] + slop, idx))

                elif format == "bedpe":
                    v.append((r['chrom'], r['start'] - slop, r['start'] + slop, r["chrom2"], r["end"] - slop,
                              r["end"] + slop, idx))

            v = ["\t".join([str(j) for j in i]) + "\n" for i in sorted(v, key=lambda x: (x[0], x[1]))]
            out.writelines(v)


def concat_dfs(cs_list):

    if isinstance(cs_list, CallSet):
        cs_list = [cs_list]
    else:
        assert isinstance(cs_list, list)

    # Add dataset info to each breaks_df
    for item in cs_list:
        item.breaks_df["caller"] = [item.caller] * len(item.breaks_df)
        item.breaks_df["dataset"] = [item.dataset] * len(item.breaks_df)
    return pd.concat([i.breaks_df for i in cs_list])


@click.command()
@click.argument('svbench_file', required=True, type=click.Path(exists=True))
@click.argument('input_file', required=False, type=click.Path())
@click.option("-o", "output", help="Output file", required=False, type=click.Path(), default="stdout",
              show_default=True)
@click.option("-c", "--col", help="Column name for model output", default="FORMAT,PROB", type=str,
              show_default=True)
@click.option("-m", "--method", help="Classifier method to use", default="predict_proba",
              type=click.Choice(["predict", "predict_proba"]), show_default=True)
def apply_model(**kwargs):
    """Apply an SVBench classifier object to SV dataset."""
    t0 = time.time()
    d = pickle.load(open(kwargs["svbench_file"], "rb"))
    d.dataset = None

    if not isinstance(d, CallSet):
        raise ValueError("SVBench file is not instance of CallSet {}".format(type(d)))
    if d.model is None:
        raise ValueError("SVBench file has no model, attribute model=None")

    new_col, add_to = None, None
    if kwargs["col"].count(",") == 1:
        add_to, new_col = kwargs["col"].split(",")
        print("Adding column {} to {}".format(new_col, add_to), file=stderr)
    else:
        raise ValueError("col argument not understood")

    if kwargs["input_file"] is None:
        kwargs["input_file"] = "-"  # assume stdin

    d.load(kwargs["input_file"])
    # Add arguments to context insert_median, insert_stdev, read_length, out_name

    if kwargs["output"] == "stdout" or kwargs["output"] == "-":
        out = stdout
    else:
        out = open(kwargs["output"], "w")

    with out:
        if kwargs["method"] == "predict_proba":
            d.predict_proba(col_name=new_col).save_as(out, new_col, add_to)

    click.echo("SVBench complete, mem={} Mb, time={} h:m:s".format(
        int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
        str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)

