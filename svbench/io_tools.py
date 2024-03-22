import ncls
from collections import defaultdict
import numpy as np
import vcf
import operator
import datetime
from importlib.metadata import version
import copy
from sys import stderr, stdin
import gzip
from io import StringIO
import networkx as nx
from collections import Counter
import pandas as pd
import os


__all__ = ["Col", "CallSet", "concat_dfs"]

svbench_version = version("svbench")

class Col:
    """This is a column parser class. The input column must be index-able by 'col' argument. Subfields are accessed
    using the 'key' argument. Values can be further encoded by providing a dict with the required mappings using the
    'encoding' argument. The 'bins' argument can be used to stratify input into predefined bins.

    :param col: The key of the primary data field
    :rtype col: str, optional
    :param key: The key of the secondary data field, optional
    :rtype key: str, optional

    """
    def __init__(self, col, key=None, encoding=None, bins=None, norm=None, op=None, thresh=None, add=None):
        self.col = col
        self.key = key
        self.encoding = encoding
        self.bins = bins
        self.norm = norm
        self.op = op
        self.thresh = thresh
        self.add = None
        if add is not None:
            self.add = add

    def __repr__(self):
        return f"svbench.Col(add={self.add} col={self.col}, key={self.key}, encoding={self.encoding}, bins={self.bins}, norm={self.norm}, op={self.op}, thresh={self.thresh})"


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
        if isinstance(a, str):
            # and a.isnumeric():
            try:
                a = float(a)
            except ValueError:
                pass
        if isinstance(b, str): # and b.isnumeric():
            try:
                b = float(b)
            except ValueError:
                pass
        try:
            v = self.opps[o](a, b)
        except TypeError:
            print(f"Failed operation using op={o}, a={a} {type(a)}, b={b} {type(b)}", file=stderr)
            quit()
        return v


class NSV:
    def __init__(self, starts, ends, ids):
        self.starts = starts
        self.ends = ends
        self.ids = ids
        self.ncls = ncls.NCLS(starts, ends, ids)

    def __reduce__(self):  # Problems pickling original ncls object, so made a new reduce method
        return self.__class__, (self.starts, self.ends, self.ids)


def get_interval_arrays(regions, slop, interval_type="breakpoint"):

    chrom_interval_start = defaultdict(list)
    chrom_interval_end = defaultdict(list)
    chrom_interval_index = defaultdict(list)

    # Make a region for each breakpoint
    if interval_type == "breakpoint":
        for c1, s, c2, e, index in regions:
            chrom_interval_start[c1].append(s - slop)
            chrom_interval_end[c1].append(s + slop)
            chrom_interval_index[c1].append(index)

            chrom_interval_start[c2].append(e - slop)
            chrom_interval_end[c2].append(e + slop)
            chrom_interval_index[c2].append(index)

    elif interval_type == "bedpe":
        for c1, s1, e1, c2, s2, e2, index in regions:
            chrom_interval_start[c1].append(s1)
            chrom_interval_end[c1].append(e1)
            chrom_interval_index[c1].append(index)

            chrom_interval_start[c2].append(s2)
            chrom_interval_end[c2].append(e2)
            chrom_interval_index[c2].append(index)

    elif interval_type == "bed":
        for c1, s, c2, e, index in regions:
            chrom_interval_start[c1].append(s)
            chrom_interval_end[c1].append(e)
            chrom_interval_index[c1].append(index)

    return ({k: np.array(v).astype(int) for k, v in chrom_interval_start.items()},
           {k: np.array(v).astype(int) for k, v in chrom_interval_end.items()},
           {k: np.array(v).astype(int) for k, v in chrom_interval_index.items()})


def make_ncls_table(chrom_interval_start, chrom_interval_end, chrom_interval_index):
    """Test doc

    :param chrom_interval_start: an array of interval start start positions
    :rtype chrom_interval_start: int, float
    :param chrom_interval_end: an array of interval start start positions
    :rtype chrom_interval_end: int, float
    :param chrom_interval_index: an array of interval start start positions
    :rtype chrom_interval_index: int, float

    :return: dict
    :rtype: dict"""

    return {k: NSV(chrom_interval_start[k], chrom_interval_end[k], chrom_interval_index[k])
            for k in chrom_interval_start}


def make_tree(regions, slop=250, interval_type="breakpoint"):
    chrom_interval_start, chrom_interval_end, chrom_interval_index = get_interval_arrays(regions, slop, interval_type)
    return make_ncls_table(chrom_interval_start, chrom_interval_end, chrom_interval_index)


def col_parser(r, col, key, first=True):
    if col == "FORMAT":
        col = r.__getattribute__("samples")
        if len(col) > 1:
            # Return the value from all the samples
            try:
                ck = [i[key] for i in col]
            except IndexError:
                raise IndexError(f"Error parsing key {key}")
            # raise ValueError("SVBench only supports one FORMAT column per vcf file, {} found", len(col))
        else:
            try:

                ck = col[0][key]
            except (IndexError, AttributeError):
                ck = None
                # raise IndexError(f"The FORMAT column is missing a name: {key}")

    else:
        try:
            col = r.__getattribute__(col)
        except KeyError:
            raise KeyError("Column argument not understood, did you mean 'INFO'?")

        if key is not None:
            if key in col:
                ck = col[key]
            else:
                ck = None

        else:
            ck = col
        if isinstance(ck, list):
            if len(ck) >= 1:  # Get first item of list or use None
                v = ck[0]
                try:
                    v = float(v)
                    return v
                except ValueError:
                    return v
            # elif len(ck) > 1:
            #     if first:
            #         ck = ck[0]
            else:
                return ck
            # else:
            #     ck = None
    return ck


def check_passed(operations, r, keep):
    passed = True
    for item in keep:
        if item.op is None:
            raise ValueError("Col.op must be set using 'keep' argument e.g. "
                             "Col('INFO', 'SU', op=eq, thresh=4)")
        ck = col_parser(r, item.col, item.key)
        if ck is None and item.col != "FILTER":
            raise KeyError("TypeError: >= 'NoneType' and 'int'. Check the svbench.Col key exists {}".format(str(keep)))

        if not operations.test(item.op, ck, item.thresh):
            passed = False
            break

    return passed


def get_strata(r, strat):

    cp = col_parser(r, strat.col, strat.key)
    if strat.add is not None:
        vals = [cp]
        current = strat.add
        while True:
            v = col_parser(r, current.col, current.key)
            if v is not None:
                vals.append(v)
            if current.add is None:
                break
        cp = sum(vals)

    if cp is None:
        cp = 0
    if isinstance(cp, list):
        cp = float(cp[0])
    return float(cp)


def parse_cols(r, item):

    encoding = item.encoding
    col = item.col
    key = item.key

    if encoding:

        if key is not None:
            p = col_parser(r, item.col, item.key)
            if p in encoding:
                item.parsed_value = encoding[p]
                return f"{col}:{key}", item.parsed_value  # encoding[p]
            elif None in encoding:
                item.parsed_value = encoding[None]
                return f"{col}:{key}", item.parsed_value
            else:

                raise ValueError("value from column={}:{} could not be found in encoding. Value was {}, of type {}".format(col, key, p, type(p)), item)

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
                item.parsed_value = encoding[p]
                return col, item.parsed_value
            else:
                item.parsed_value = encoding[None]
                return col, item.parsed_value

    else:
        if key is not None:
            p = col_parser(r, col, key, first=False)
            if isinstance(p, list):
                d = dict()
                for index, v in enumerate(p):
                    d[f"{col}:{key}:{index}"] = v
                return d
            else:
                item.parsed_value = p
                return f"{col}:{key}", item.parsed_value
        else:
            try:
                item.parsed_value = r[col]
            except TypeError:
                item.parsed_value = r.__getattribute__(col)
            return col, item.parsed_value


def parse_cols_list(r, parse_list, data):
    for item in parse_list:
        d = parse_cols(r, item)
        if isinstance(d, dict):
            data.update(d)
        else:
            data[d[0]] = d[1]
    return data


def parse_all_cols(r):
    all_cols = []
    for k in r.INFO.keys():
        all_cols.append(Col("INFO", k))
    for k in r.FORMAT.split(":"):
        all_cols.append(Col("FORMAT", k))
    return parse_cols_list(r, all_cols, {})


def check_args(stratify, weight_field, keep, other_cols, size_range):

    if size_range is not None:
        if size_range[0] is not None and size_range[1] is not None:
            assert size_range[0] < size_range[1]

    for item in (keep, weight_field, stratify):
        if item is None:
            continue
        if isinstance(item, list):
            for item2 in item:
                assert isinstance(item2, Col)
        else:
            assert isinstance(item, Col)

    if isinstance(other_cols, list):
        for item2 in other_cols:
            assert isinstance(item2, Col)
    # elif other_cols != "all":
    #     assert isinstance(other_cols, Col)

    if stratify is not None and (stratify.bins is None or not hasattr(stratify.bins, '__iter__')):
        raise ValueError(
            "Stratify must have an iterable for 'bins' argument e.g. Col('FORMAT', 'DV', bins=range(0, 60, 5))")


class CallSet:
    """
    This class holds instructions for parsing an input file containing variants. Supported formats include
    'vcf', 'bed', 'bedpe', 'csv'.
    Raw data may optionally be saved within the object.

    :param dataset: A name for the input dataset (optional)
    :type dataset: str
    :param caller: The name of the variant caller used (optional)
    :type caller: str
    :param kwargs: Allows CallSet attributes to be set during initialization - if any kwargs are not found in the \
    self.default_params.keys() list, then the key: value pair is added to the self.required dictionary. The purpose of \
    this feature is to allow the setting of required arguments, which can be useful when creating new class instances.
    :returns: CallSet instance
    """
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
        self.weight_field = None
        self.no_translocations = True
        self.allowed_svtypes = None
        self.keep = None
        self.stratify = None
        self.stratify_range = None
        self.allowed_chroms = None
        self.min_size = None
        self.max_size = None
        self.soft_size_filter = None
        self.other_cols = None
        self.model = None
        self.slop = None
        self.break_cols = "chrA,posA,chrB,posB"
        self.sep = ","
        self.drop_first_row = False
        self.svtype_field = "svtype"
        self.scores = None
        self.size_scores = None
        self.false_negative_indexes = None
        self.kind = None
        self.new_col = None
        self.add_to = None
        self.add_chr_prefix = True
        self.id_field = None
        self.temp = None
        self.include_bed = None
        self.include_if = "both"
        self.load_genotype = True
        self.gt_scores = None

        # Track defaults for persistence
        self.default_params = {k: v for k, v in self.__dict__.items() if k in {"bedpe",
                                                                               "weight_field",
                                                                               "no_translocations",
                                                                               "allowed_svtypes",
                                                                               "keep",
                                                                               "stratify",
                                                                               "stratify_range",
                                                                               "allowed_chroms",
                                                                               "min_size",
                                                                               "max_size",
                                                                               "soft_size_filter",
                                                                               "other_cols",
                                                                               "slop",
                                                                               "break_cols",
                                                                               "sep",
                                                                               "drop_first_row",
                                                                               "svtype_field"}}

        self.meta_data = {"svbench": svbench_version,
                          "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

        # Update any provided kwargs
        if kwargs:
            self.set_properties(kwargs)

    def __add__(self, other):
        if not isinstance(other, CallSet):
            raise ValueError("Calling add on {}, should be instance of CallSet".format(type(other)))

        n = self.copy()

        n.breaks_df["source"] = [n.caller] * len(n.breaks_df)
        other_df = other.breaks_df.copy(deep=True)
        other_df["source"] = [other.caller] * len(other_df)
        n.breaks_df = pd.concat([n.breaks_df, other_df])
        n.breaks_df.reset_index(inplace=True)
        n.tree = None
        return n

    def __sub__(self, other):
        if not isinstance(other, CallSet):
            raise ValueError("Calling add on {}, should be instance of CallSet".format(type(other)))
        if self.tree is None or other.tree is None:
            raise ValueError("Call add_intervals on both CallSet's before subtracting")

        matching_indexes = quantify(other, self, good_indexes_only=True)
        n = self.copy()
        n.breaks_df = n.breaks_df[~np.array(matching_indexes)]
        n.tree = None
        return n

    def __len__(self):
        return len(self.breaks_df) if self.breaks_df is not None else 0

    def intersection(self, other):
        return self.__iand__(other)

    def __iand__(self, other):
        if not isinstance(other, CallSet):
            raise ValueError("Calling add on {}, should be instance of CallSet".format(type(other)))
        if self.tree is None or other.tree is None:
            raise ValueError("Call add_intervals on both CallSet's before subtracting")

        matching_indexes = quantify(other, self, good_indexes_only=True)
        n = self.copy()
        n.breaks_df = n.breaks_df[np.array(matching_indexes)]
        n.tree = None
        return n

    def add_intervals(self, slop=250, interval_type="breakpoint"):
        """Adds intervals to loaded data defined in self.breaks_df.

        :param slop: The distance to add and subtract from each variant position.
        :type slop: int
        :returns: self
        :rtype: svbench.CallSet"""
        if interval_type not in ("breakpoint", "bed", "bedpe"):
            raise ValueError("interval_type must be one of: 'breakpoint', 'bed', 'bedpe'")

        self.slop = slop
        df = self.breaks_df
        if df is None or len(df) == 0:
            raise ValueError("breaks_df not found")
        self.tree = make_tree(list(zip(df["chrom"], df["start"], df["chrom2"], df["end"], df.index)), slop=slop,
                              interval_type=interval_type)
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

    def copy(self):
        """Create a new deep-copy of this object

        :returns: CallSet instance
        :rtype: svbench.CallSet"""
        return copy.deepcopy(self)

    def new(self, **kwargs):
        n = self.copy()
        n.set_properties(kwargs)
        return n

    def query(self, query_expression, inplace=False, engine="python"):
        if not inplace:
            cs = self.copy()
        else:
            cs = self
        l_before = len(cs.breaks_df)
        cs.breaks_df = cs.breaks_df.query(query_expression, engine=engine)
        print("Filtered by expression, caller={}, dataset={} rows before {}, after {}".format(self.caller, self.dataset,
                                                                                          l_before,
                                                                                          len(cs.breaks_df)), file=stderr)
        return cs

    def filter_by_size(self, min_size=None, max_size=None, inplace=False, soft=False, keep_translocations=False):
        """Filter the loaded data defined in self.breaks_df by size (base-pairs).

        :param min_size: The minimum size threshold for the variant, size < min_size
        :type min_size: int
        :param max_size: The maximum size threshold for the variant, size > max_size
        :type max_size: int
        :param inplace: If set to True then filtering is performed inplace
        :type inplace: bool
        :keep_translocations: False means translocations calls will be dropped
        :return: CallSet instance
        :rtype: svbench.CallSet
        """
        if inplace:
            cs = self
        else:
            cs = self.copy()
        done = True
        if cs.min_size != min_size:
            cs.min_size = min_size
            done = False
        if cs.max_size != max_size:
            cs.max_size = max_size
            done = False
        if done or (min_size is None and max_size is None):  # Nothing to be done
            return cs

        if min_size is None:
            min_s = -1
        else:
            min_s = min_size
        if max_size is None:
            max_s = 1e12
        else:
            max_s = max_size

        l_before = len(cs.breaks_df)
        df = cs.breaks_df

        if keep_translocations:
            size_filter = np.array([True if ((svlen >= min_s and svlen < max_s) or chrom1 != chrom2 or svlen is None or svlen != svlen) else False for svlen, chrom1, chrom2 in zip(df['svlen'], df['chrom'], df['chrom2'])])
        else:
            size_filter = np.array([True if ((svlen >= min_s and svlen < max_s) or svlen is None or svlen != svlen) else False for svlen, chrom1, chrom2 in zip(df['svlen'], df['chrom'], df['chrom2'])])

        df["size_filter_pass"] = size_filter
        if not soft:
            df = df[size_filter]
        print("Filtered by min_size={}, max_size={}, caller={}, dataset={} rows before {}, after {}".format(
            min_size, max_size, self.caller, self.dataset,l_before, size_filter.sum()),
              file=stderr)

        cs.breaks_df = df
        return cs

    def filter_by_svtype(self, svtype_set, inplace=False):
        l_before = len(self.breaks_df)
        bad_i = set([])
        if isinstance(svtype_set, str):
            svtype_set = set(tuple(svtype_set.split(",")))
        for idx, svtype in zip(self.breaks_df.index, self.breaks_df.svtype):
            if svtype not in svtype_set:
                bad_i.add(idx)
        print("Filtered by svtype, caller={}, dataset={} rows before {}, after {}".format(self.caller, self.dataset,
                                                                                          l_before,
                                                                                          l_before - len(bad_i)),
              file=stderr)
        if not inplace:
            s = self.copy()
            s.breaks_df = s.breaks_df.drop(bad_i)
            return s
        self.breaks_df = self.breaks_df.drop(bad_i)

        return self

    def filter_include_bed(self, include_bed_path, inplace=False):
        if not inplace:
            v = self.copy()
        else:
            v = self
        l_before = len(v.breaks_df)
        bad_i = set([])
        if include_bed_path:
            include_bed = CallSet().load_bed(path=include_bed_path, bedpe=False).add_intervals(interval_type="bed")
            ol_tree = include_bed.tree
            df = v.breaks_df
            if df is not None and ol_tree is not None:
                for index, chrom, start, chrom2, end in zip(df.index, df.chrom, df.start, df.chrom2, df.end):
                    if ol_tree:
                        if chrom not in ol_tree:
                            bad_i.add(index)
                        else:
                            if not any(ol_tree[chrom].ncls.find_overlap(start, start + 1)) or not any(ol_tree[chrom2].ncls.find_overlap(end, end + 1)):
                                bad_i.add(index)
            v.breaks_df = v.breaks_df.drop(bad_i)

        print("Filtered by include_bed, caller={}, dataset={} rows before {}, after {}".format(v.caller, v.dataset,
                                                                                          l_before,
                                                                                          l_before - len(bad_i)),
              file=stderr)
        return v

    def score_table(self):
        if self.scores is not None:
            print(f"Score table caller={self.caller} against dataset={self.dataset}", file=stderr)
            print(pd.DataFrame.from_records([self.scores], index=None).to_string(), file=stderr)

    def set_strata(self, stratify_col, stratify_range):

        if not hasattr(stratify_range, '__iter__'):
            raise ValueError("stratify_range must be an iterable")

        self.breaks_df["strata"] = self.breaks_df[stratify_col]
        self.stratify_range = stratify_range
        return self

    def load(self, path):
        """Load variants from the file path. For this function to work the CallSet instance must have the self.kind \
        attribute set.

        :param path: File path for input data
        :type path: str
        :return: CallSet instance
        :rtype: svbench.CallSet
        """
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

    def parse_str_cols(self, cols):
        if type(cols) == list and all(type(i) == str for i in cols):
            c = []
            for k in cols:
                if ":" in k:
                    c.append(Col(*k.split(":")))
                else:
                    c.append(Col(k))
            return c
        else:
            return cols

    def load_vcf(self, path, weight_field=None,
                 allowed_svtypes=None, keep=None, stratify=None, allowed_chroms=None,
                 min_size=None, max_size=None, soft_size_filter=False,
                 other_cols=None, include_bed=None, include_if="both", load_genotype=True):
        """Load variants from the vcf file path.

        :param path: The path to the vcf input file
        :type path: str
        :param weight_field: The field used to weight a variant, useful for breaking ties between similar variants when \
        benchmarking
        :type weight_field: svbench.Col
        :type no_translocations: bool
        :param allowed_svtypes: The types of SVs allowed, comma seperated string e.g. "DEL", or "DEL,DUP,INV" etc
        :type allowed_svtypes: str
        :param keep: A list of filtering operations to perform on input data. e.g. keep=[Col("INFO", "SU", "ge" 3)] \
        would keep variants if the value found in INFO --> SU was greater than or equal to 3
        :type keep: list
        :param stratify:
        :param allowed_chroms:
        :param min_size:
        :param max_size:
        :param soft_size_filter:
        :param other_cols:
        :param include_bed:
        :param include_if:
        :return:
        """

        self.path = path
        self.kind = "vcf"

        other_cols = self.parse_str_cols(other_cols)

        check_args(stratify, weight_field, keep, other_cols, (min_size, max_size))

        self.set_args(locals())  # Overwrite default params

        ol_tree = None
        if include_bed is not None:
            if isinstance(include_bed, CallSet):
                assert include_bed.tree is not None
                assert include_if in ("both", "single")
                ol_tree = include_bed.tree
            elif isinstance(include_bed, str):
                v = CallSet(dataset=self.dataset)
                include_bed = v.load_bed(path=include_bed, bedpe=False).add_intervals()
                ol_tree = include_bed.tree
            else:
                raise ValueError("include_bed must be of type svbench.CallSet or PATH")

        # Load instance parameters, these may have been set previously
        add_weight = self.weight_field is not None
        weight_field = self.weight_field
        no_translocations = self.no_translocations
        allowed_svtypes = self.allowed_svtypes
        keep = self.keep
        stratify = self.stratify
        allowed_chroms = self.allowed_chroms
        min_size = self.min_size
        max_size = self.max_size
        soft_size_filter = self.soft_size_filter
        other_cols = self.other_cols

        operations = Operate()

        res = []
        if allowed_svtypes is not None:
            allowed_svtypes = set(allowed_svtypes.split(","))

        if allowed_chroms:
            if not isinstance(allowed_chroms, set):
                allowed_chroms = set(allowed_chroms)

        new_cols = []
        unique_ids = set([])
        ignore_mates = set([])

        if path in "-stdin":
            temp = StringIO()
            for line in stdin:
                temp.write(line)
            temp.seek(0)
            self.temp = temp  # Used when writing output

            reader = vcf.Reader(fsock=temp)
        else:
            assert os.path.exists(path)
            reader = vcf.Reader(filename=path)

        if other_cols == "all":
            l = []
            for k in reader.infos.keys():
                l.append(Col("INFO", k))
            for k in reader.formats.keys():
                l.append(Col("FORMAT", k))
            other_cols = l

        #
        vcf_index = 0
        not_in_include = 0
        while True:
            try:
                r = next(reader)
                vcf_index += 1
            except StopIteration:
                break
            except IndexError:
                raise IndexError("last record was", str(r))

            ol_start = False
            ol_end = False

            chrom = r.CHROM
            if isinstance(chrom, int) or (isinstance(chrom, str) and chrom[0] != "c"):
                chrom = "chr" + str(chrom)

            if allowed_chroms is not None and chrom not in allowed_chroms:
                continue

            if "#" in chrom:
                continue

            if "MATEID" in r.INFO:
                mate = r.INFO["MATEID"]
                if isinstance(mate, str):
                    if mate in ignore_mates:
                        continue
                    ignore_mates.add(mate)
                elif isinstance(mate, list):
                    skip_this = False
                    for m in mate:
                        if m in ignore_mates:
                            skip_this = True
                            break
                        ignore_mates.add(m)
                    if skip_this:
                            continue

            start = int(r.POS)

            if ol_tree and chrom in ol_tree:
                if any(ol_tree[chrom].ncls.find_overlap(start, start + 1)):
                    ol_start = True
                if not ol_start and include_if == "both":
                    not_in_include += 1
                    continue

            if "SVTYPE" in r.INFO:
                svtype = r.INFO["SVTYPE"] if not isinstance(r.INFO["SVTYPE"], list) else r.INFO["SVTYPE"][0]
            else:
                if isinstance(r.REF, list):
                    lr = len(r.REF[0])
                else:
                    lr = len(r.REF)
                if isinstance(r.ALT, list):
                    la = len(r.ALT[0])
                else:
                    la = len(r.ALT)
                if lr > 1 or la > 1:
                    if lr > la:
                        r.INFO["SVTYPE"] = "DEL"
                        svtype = "DEL"
                    else:
                        r.INFO["SVTYPE"] = "INS"
                        svtype = "INS"
                else:
                    r.INFO["SVTYPE"] = "BND"
                    svtype = "BND"

            try:
                if allowed_svtypes is not None and svtype not in allowed_svtypes:
                    continue
            except TypeError:
                if isinstance(svtype, list):
                    svtype = svtype[0]
                    if allowed_svtypes is not None and svtype not in allowed_svtypes:
                        continue

            if keep is not None and not check_passed(operations, r, keep):
                continue

            svlen = -1
            if isinstance(r.ID, str) and r.ID.endswith("_2"):  # Skip second part of BND, or multirow record
                continue
            if svtype == "BND":
                done = False
                if "CHR2" in r.INFO:
                    chrom2 = r.INFO["CHR2"]
                    if "CHR2_POS" in r.INFO:
                        end = r.INFO["CHR2_POS"]
                        done = True
                    if "END" in r.INFO:
                        if isinstance(r.INFO["END"], str):
                            try:
                                if "[" in r.ALT[0]:
                                    chrom2, pos2 = [i for i in r.ALT[0].split("[") if ":" in i][0].split(":")
                                    end = int(pos2)
                                    done = True
                                elif "]" in r.ALT[0]:
                                    chrom2, pos2 = [i for i in r.ALT[0].split("]") if ":" in i][0].split(":")
                                    end = int(pos2)
                                    done = True
                            except:
                                pass
                        elif isinstance(r.INFO["END"], int):
                            end = r.INFO["END"]
                            done = True
                if not done:
                    if r.ALT[0] is None:
                        continue
                    try:
                        chrom2 = r.ALT[0].chr
                        end = r.ALT[0].pos
                        done = True
                    except:
                        pass
                    if not done:
                        try:
                            chrom2 = r.CHROM
                            end = r.end
                        except AttributeError:
                            print("AttributeError parsing", r.ID, file=stderr)
                            continue
                        if end is None:
                            end = start + 1
                        else:
                            end = r.POS
                        if chrom2 is None:
                            chrom2 = chrom  # give up
                if chrom.startswith("chr") and not chrom2.startswith("chr"):
                    chrom2 = "chr" + chrom2
            else:
                chrom2 = chrom
                if "CHR2" in r.INFO:
                    chrom2 = r.INFO["CHR2"]
                    if chrom2[0] != "c":
                        chrom2 = "chr" + chrom2

                done = False
                if "END" in r.INFO or "CHR2_POS" in r.INFO:
                    if "CHR2_POS" in r.INFO:
                        end = r.INFO["CHR2_POS"]
                        if isinstance(end, list):
                            end = end[0]
                        end = int(end)
                        done = True
                    else:
                        end = r.INFO["END"]
                        if isinstance(end, list):
                            end = end[0]
                        end = int(end)
                        done = True
                elif svtype in ("DEL", "DUP", "INV") and "SVLEN" in r.INFO:
                    svlen = r.INFO["SVLEN"]
                    if isinstance(svlen, list):
                        svlen = svlen[0]
                    end = start + svlen
                    done = True
                else:  # Try and use ALT / REF lengths
                    if r.INFO["SVTYPE"] == "DEL" or (isinstance(r.INFO["SVTYPE"], list) and r.INFO["SVTYPE"][0] == "DEL"):
                        svlen = len(r.REF) if r.ALT is not isinstance(r.REF, list) else r.REF[0]
                        end = r.POS + svlen
                        done = True
                    elif r.INFO["SVTYPE"] == "INS" or (isinstance(r.INFO["SVTYPE"], list) and r.INFO["SVTYPE"][0] == "INS"):
                        svlen = len(r.ALT) if r.ALT is not isinstance(r.ALT, list) else r.ALT[0]
                        end = r.POS + svlen
                        done = True

                if not done:
                    print("Warning: could not parse record", r, r.INFO, file=stderr)
                    continue

            if no_translocations and chrom != chrom2:
                continue

            if ol_tree and chrom2 in ol_tree:
                if any(ol_tree[chrom2].ncls.find_overlap(end, end + 1)):
                    ol_end = True
                if not ol_end and include_if == "both":
                    not_in_include += 1
                    continue
                elif not ol_start and not ol_end and include_if == "single":
                    not_in_include += 1
                    continue

            size_filter = True

            if chrom == chrom2:
                if "SVLEN" in r.INFO:
                    svlen = r.INFO["SVLEN"]
                    if isinstance(svlen, list):
                        svlen = svlen[0]

                    if svlen is None:
                        svlen = None
                    else:
                        svlen = abs(svlen)

                else:
                    svlen = abs(end - start)
                    if svlen < 2 and "SVTYPE" in r.INFO and r.INFO["SVTYPE"] == "INS":
                        svlen = min_size

                if min_size is not None and svlen < min_size:
                    if not soft_size_filter:
                        continue
                    else:
                        size_filter = False
                if max_size is not None and svlen >= max_size:
                    if not soft_size_filter:
                        continue
                    else:
                        size_filter = False

            if allowed_chroms is not None and chrom2 not in allowed_chroms:
                continue

            w = None
            if add_weight:
                w = parse_cols(r, weight_field)
                if isinstance(w, dict):
                    raise ValueError("Multiple values in weight field, choose a field with a single value")
                wcol, w = w

            if r.ID in unique_ids or r.ID is None:
                r_id = f"{r.CHROM}:{r.POS}.line_{vcf_index}"
                unique_ids.add(r_id)
            else:
                r_id = r.ID
                unique_ids.add(r_id)

            d = {"end": end, "start": start, "chrom": chrom, "chrom2": chrom2, "w": w, "id": r_id, "svtype": svtype,
                 "size_filter_pass": size_filter, "svlen": abs(svlen) if svlen else svlen, "filter": r.FILTER}

            if stratify is not None:
                d["strata"] = get_strata(r, stratify)

            if other_cols:
                parsed = parse_cols_list(r, other_cols, {})
                if not new_cols:
                    new_cols = list(parsed.keys())
                d.update(parsed)

            samps = r.__getattribute__("samples")
            if load_genotype and len(samps) == 1:
                # if len(samps) > 1:
                    # raise ValueError("Cannot parse genotype for multi-sample vcf, set load_genotype=False")
                done = False
                try:
                    d["GT"] = str(samps[0]["GT"])
                    done = True
                except IndexError:
                    pass

                if not done:
                    # try:
                    #     d["GT"] = r.INFO["GT"]
                    # except KeyError:
                    #     pass
                    d["GT"] = "NA"
            else:
                d["GT"] = "NA"

            res.append(d)

        if len(res) == 0:
            print(f"Warning: empty vcf {path}", file=stderr)
            return self

        df = pd.DataFrame.from_records(res)

        # Normalize new columns, might be more than one column added per Col
        if other_cols is not None and other_cols != "all":
            item_key = {(k.col, k.key): k for k in other_cols if k.norm is not None}
            for k, v in item_key.items():
                for ec in new_cols:
                    ecs = ec.split(":")
                    if ecs[0] == k or (ecs[0], ecs[1]) == k:
                        df[ec] = v.norm(df[ec], self.kwargs)

        # Order df
        base_cols = ["chrom", "start", "chrom2", "end", "svtype", "w", "strata", "id", "size_filter_pass", "svlen", "filter"]
        if load_genotype:
            base_cols.append("GT")
            df["GT"] = [i.replace("|", "/") if isinstance(i, str) else i for i in df["GT"]] #df["GT"].str.replace("|", "/")
        df = df[[i for i in base_cols if i in df.columns] + new_cols]

        self.breaks_df = df
        self.extra_cols = new_cols

        if stratify is not None:
            self.stratify_range = stratify.bins

        print(f"dataset={self.dataset}, caller={self.caller}, loaded rows: {len(df)}", file=stderr)
        print("columns:", list(self.breaks_df.columns), file=stderr)
        print("n not in include: {}".format(not_in_include), file=stderr)
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
            if np.issubdtype(df['chrom'].dtype, np.number) or str(df['chrom'].iloc[0]).isdigit() or df["chrom"].iloc[0][0] != "c":
                df["chrom"] = ["chr" + str(i) for i in df["chrom"]]
            df["chrom2"] = df["chrom"]
            df["start"] = df["start"].astype(int)
            df["end"] = df["end"].astype(int)

        df["strata"] = [None] * len(df)
        df["id"] = df.index
        df["size_filter_pass"] = [True] * len(df)
        df["w"] = [1] * len(df)
        svlen = []
        for chrom1, chrom2, start, end in zip(df["chrom"], df["chrom2"], df["start"], df["end"]):
            if chrom1 == chrom2:
                svlen.append(abs(end - start))
            else:
                svlen.append(-1)
        df["svlen"] = svlen

        self.breaks_df = df
        print(f"dataset={self.dataset}, caller={self.caller}, loaded rows: {len(df)}", file=stderr)

        return self

    def load_csv(self, path, break_cols="chrA,posA,chrB,posB", sep=",", weight_field=None,
                 allowed_svtypes=None, keep=None, svtype_field="svtype", no_translocations=True, id_field=None,
                 allowed_chroms=None, stratify=None,
                 min_size=None, max_size=None, soft_size_filter=False,
                 other_cols=None, drop_first_row=False,
                 add_chr_prefix=True):

        self.kind = "csv"

        check_args(stratify, weight_field, keep, other_cols, (min_size, max_size))

        self.set_args(locals())

        # Load instance parameters, these may have been set previously
        break_cols = self.break_cols
        sep = self.sep
        add_weight = self.weight_field is not None
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

        if isinstance(path, str):
            df_in = pd.read_csv(path, sep=sep, index_col=None, comment="#")
        else:
            # Assume pandas dataframe
            if not isinstance(path, pd.DataFrame):
                raise ValueError("Only pandas DataFrame objects of 'str' arguments can be accepted")
            df_in = path

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

            if n[:3] == "chr" and add_chr_prefix and df_in[col].iloc[0][:3] != "chr":
                df[n] = df_in[col].astype(str).add_prefix("chr")
            else:
                df[n] = df_in[col]

        if id_field is not None:
            df["id"] = df_in[id_field]
        else:
            df["id"] = df.index

        if allowed_chroms:
            if not isinstance(allowed_chroms, set):
                allowed_chroms = set(allowed_chroms)
            df = df[df["chrom"].isin(allowed_chroms) & df["chrom2"].isin(allowed_chroms)]

        if no_translocations:
            df = df[df["chrom"] == df["chrom2"]]

        df["size_filter_pass"] = [True] * len(df)
        if min_size is not None or max_size is not None:
            drop = []
            for idx, start, end, chrom1, chrom2 in zip(df.index, df["start"], df["end"], df["chrom"], df["chrom2"]):
                if chrom1 == chrom2:
                    s = abs(end - start)
                    if min_size is not None and s < min_size:
                        drop.append(idx)
                    elif max_size is not None and s >= max_size:
                        drop.append(idx)

            if not soft_size_filter:
                df.drop(index=drop, inplace=True)
            else:
                df["size_filter_pass"] = [i not in drop for i in df.index]

        if stratify is not None:
            df["strata"] = df_in[stratify.col].loc[df.index]
            self.stratify_range = stratify.bins

        if add_weight:
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

    def write_to_vcf(self, path, new_col=None, add_to="FORMAT", drop_missing=False):

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

        if new_col is not None:
            col_vals = {idx: v for idx, v in zip(self.breaks_df.id, self.breaks_df[new_col])}
        else:
            col_vals = {idx: None for idx in self.breaks_df.id}

        replace_val = False
        write_h_line = True
        write_h = True
        header_lines = []
        last_add_index = 0
        position = 0
        print("Writing vcf to", path, file=stderr)
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
                    r_id = l[2]

                    if r_id in col_vals:
                        vf = col_vals[r_id]
                        if vf is None:
                            val = None
                        elif vf == 1:
                            val = "1"
                        elif vf == 0:
                            val = "0"
                        else:
                            val = f"{vf:.4f}"

                    elif f"{l[0]}:f{l[1]}" in col_vals:
                        vf = col_vals[f"{l[0]}:f{l[1]}"]
                        if vf is None:
                            val = None
                        elif vf == 1:
                            val = "1"
                        elif vf == 0:
                            val = "0"
                        else:
                            val = f"{vf:.4f}"
                    elif drop_missing:
                        continue
                    else:
                        val = "0"
                    if new_col is None:
                        f_out.writelines("\t".join(l))
                    elif add_to == "FORMAT":
                        vidx = None
                        if replace_val:
                            l8s = l[8].split(":")
                            vidx = [idx for idx, vv in enumerate(l8s) if vv == new_col][0]
                            # if vidx == len(l8s)-1:
                            #     vidx = None  # Place at end
                        else:
                            if new_col is not None:
                                l[8] += ":" + new_col
                        # Value is added to each record in the row

                        for i in range(9, len(l)):
                            # if i == len(l) - 1 and vidx is None:
                            #     l[i] = l[i].strip() + ":" + val + "\n"
                            if vidx is None:
                                 l[i] = l[i] + ":" + val
                            else:
                                lis = l[i].split(":")
                                lis[vidx] = val
                                l[i] = ":".join(lis)
                                if i == len(l) - 1:
                                    l[i] += "\n"
                        f_out.writelines("\t".join(l))
                    elif add_to == "INFO":
                        l[7] += f";{new_col}={val}"
                        f_out.writelines("\t".join(l))
                    else:
                        raise ValueError("Add to INFO/FORMAT only")

                    position += 1

    def save_as(self, path, new_col=None, add_to=None, format_t=None):
        if new_col is None and self.new_col is not None:
            new_col = self.new_col
        if add_to is None and self.add_to is not None:
            add_to = self.add_to

        if format_t is None and self.kind == "vcf":
            self.write_to_vcf(path, new_col, add_to)
        elif format_t is None and self.kind == "csv":
            raise ValueError("Not implemented currently")

        else:
            raise ValueError(self.kind, "not in {csv, vcf}")

    @staticmethod
    def _sort_func(x):
        return x[0], x[1]

    def save_intervals(self, path, format_t="bed"):
        if format_t not in {"bed", "bedpe"}:
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
                if format_t == "bed":

                    v.append((r['chrom'], r['start'] - slop, r['start'] + slop, idx, r["id"]))
                    v.append((r['chrom2'], r['end'] - slop, r['end'] + slop, idx, r["id"]))

                elif format_t == "bedpe":
                    v.append((r['chrom'], r['start'] - slop, r['start'] + slop, r["chrom2"], r["end"] - slop,
                              r["end"] + slop, idx,  r["id"]))

            v = ["\t".join([str(j) for j in i]) + "\n" for i in v]
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


def quantify(ref_data, data, force_intersection=False, reciprocal_overlap=0., show_table=True, stratify=False,
             good_indexes_only=False, ref_size_bins=(30, 50, 500, 5000, 260000000), allow_duplicate_tp=True,
             pct_size=0.05, ignore_svtype=True, min_ref_size=20, max_ref_size=None, dups_and_ins_equivalent=False):

    # Build a Maximum Bipartite Matching graph:
    # https://www.geeksforgeeks.org/maximum-bipartite-matching/

    G = nx.Graph()

    ref_bedpe = ref_data.breaks_df
    dta = data.breaks_df

    tree = ref_data.tree
    if tree is None:
        raise ValueError("Interval tree has not been created, call add_intervals first")

    # # Link query calls to reference calls
    if dta is None:
        ts = [{"Total": None, "Ref": len(ref_bedpe), "DTP": None, "TP": None, "FP": None, "FN": None, "T >=": None,
               "Precision": None, "Recall": None, "F1": None, "Duplication": None, "quantified": None, "strata": None}]
        # return data
        data.scores = pd.DataFrame.from_records(ts)[["T >=", "Ref", "Total", "TP", "FP", "DTP", "FN", "Duplication", "Precision",
                                                     "Recall", "F1"]]
        data.false_negative_indexes = ref_bedpe.index
        return

    for query_idx, chrom, start, chrom2, end, svtype, w, d_id in zip(dta.index, dta["chrom"], dta["start"], dta["chrom2"], dta["end"],
                                                       dta["svtype"], dta["w"], dta["id"]):
        if chrom == chrom2 and start == end:
            end += 1  # prevent 0 width interval
        chrom, start, chrom2, end = sv_key(chrom, start, chrom2, end)

        ol_start = intersecter(tree, chrom, start, start + 1)
        # if d_id == '54109':
        #     print('start', ol_start, (start, start+1), file=stderr)
        #     quit()
        if not ol_start:
            # print("FALSE1", chrom, start, chrom2, end, svtype, file=stderr)
            continue

        ol_end = intersecter(tree, chrom2, end, end + 1)
        # if d_id == '54109':
        #     print('end', ol_end, file=stderr)
        #     quit()
        if not ol_end:
            # print("FALSE2", chrom, start, chrom2, end, svtype, file=stderr)
            continue

        # if d_id == '54109':
        #     print(ol_start, file=stderr)
        #     print(ol_end, file=stderr)
        #     quit()

        # Get the ref_data index
        common_idxs = set([i[2] for i in ol_start]).intersection([i[2] for i in ol_end])
        if len(common_idxs) == 0:
            continue

        # Choose an index by highest weight/lowest total distance, meeting reciprocal_overlap threshold
        min_d = 1e12
        chosen_index = None
        ref_chrom, ref_start, ref_chrom2, ref_end = None, None, None, None
        for index in common_idxs:
            try:
                ref_row = ref_bedpe.loc[index]
            except IndexError:
                raise ("Index error", index, " Try re-setting intervals")
            ref_chrom, ref_start, ref_chrom2, ref_end = sv_key(ref_row["chrom"], ref_row["start"],
                                                               ref_row["chrom2"], ref_row["end"])
            # Make sure chromosomes match
            if chrom != ref_chrom or chrom2 != ref_chrom2:
                continue
            if not ignore_svtype and ref_row["svtype"] != svtype:
                if dups_and_ins_equivalent and svtype in ("DUP", "INS") and ref_row["svtype"] in ("DUP", "INS"):
                    pass
                if svtype == "BND" or ref_row["svtype"] == "BND":  # lenient match
                    pass
                else:
                    continue
            # If intra-chromosomal, check reciprocal overlap
            if chrom == chrom2:
                if "svlen" in ref_row:
                    ref_size = ref_row["svlen"] + 1e-6
                else:
                    ref_size = ref_end - ref_start + 1e-3
                if min_ref_size is not None:
                    if ref_size < min_ref_size:
                        continue
                if max_ref_size is not None:
                    if ref_size >= max_ref_size:
                        continue
                if "svlen" in dta:
                    query_size = dta["svlen"].loc[query_idx] + 1e-6
                else:
                    query_size = end - start + 1e-6
                ol = float(max(0, min(end, ref_end) - max(start, ref_start)))
                if pct_size > 0:
                    pct = min(ref_size, query_size) / max(ref_size, query_size)
                    if pct < pct_size:
                        continue
                if reciprocal_overlap > 0:
                    if (ol / query_size < reciprocal_overlap) and (ol / ref_size < reciprocal_overlap):
                        continue
                if force_intersection and ol == 0:
                    continue
            dis = abs(ref_start - start) + abs(ref_end - end)
            if dis < min_d:
                min_d = dis
                chosen_index = index

        if chosen_index is not None:
            G.add_edge(('t', chosen_index), ('q', query_idx), dis=min_d, weight=w)
            # print(chosen_index, query_idx, (ref_chrom, ref_start), (ref_chrom2, ref_end), chrom, start, chrom2, end, svtype, file=stderr)

        # if start == 57861455:
        #     print(common_idxs)
        #     print(chosen_index)
        #     print(chrom, ref_chrom, chrom2, ref_chrom2, dis < min_d, index)
        #     quit()

    good_idxs = {}
    duplicate_idxs = {}
    # Make partitions bipartite-matching i.e. one reference call matched to one query call
    for sub in nx.connected_components(G):
        sub = list(sub)
        if len(sub) == 2:  # One ref matches one query, easy case
            if sub[0][0] == "q":
                good_idxs[sub[0][1]] = sub[1][1]
            else:  # t first
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
                print("Maximum bipartite matching not implemented", file=stderr)
                quit()
        continue
    df = data.breaks_df
    index = df.index  # Use the original IDs?
    if good_indexes_only:
        return [i in good_idxs for i in index]

    missing_ref_indexes = set(ref_bedpe.index).difference(set(good_idxs.values()))
    if allow_duplicate_tp:
        duplicate_idxs = {k: v for k, v in duplicate_idxs.items() if k not in good_idxs}
    else:
        duplicate_idxs = {}

    df["TP"] = [i in good_idxs for i in index]
    df["DTP"] = [i in duplicate_idxs for i in index]
    df["ref_index"] = [good_idxs[i] if i in good_idxs else duplicate_idxs[i] if i in duplicate_idxs
                                  else None for i in index]
    df["FP"] = [not i and not j for i, j in zip(df["DTP"], df["TP"])]
    if "svlen" in ref_bedpe:
        df["ref_size"] = [abs(ref_bedpe["svlen"].loc[good_idxs[i]]) if
                          (i in good_idxs and dta.loc[i]["chrom"] == dta.loc[i]["chrom2"]) else None
                          for i in index]

    else:
        df["ref_size"] = [abs(ref_bedpe["end"].loc[good_idxs[i]] - ref_bedpe["start"].loc[good_idxs[i]]) if
                                      (i in good_idxs and dta.loc[i]["chrom"] == dta.loc[i]["chrom2"]) else None
                                      for i in index]

    n_in_ref = len(ref_bedpe)
    if stratify and data.stratify_range is not None:
        rng = data.stratify_range
    else:
        rng = [None]

    data.breaks_df = df
    dta = df
    ts = []
    for threshold in rng:
        if threshold is None or len(dta) == 0:
            t = {"Total": len(dta),
                 "Ref": n_in_ref,
                 "DTP": np.sum(np.in1d(dta["DTP"], True)),
                 "TP": np.sum(np.in1d(dta["TP"], True)),
                 "FP": np.sum(np.in1d(dta["FP"], True)),
                 "FN": len(missing_ref_indexes),
                 "T >=": None,
                 "Caller": data.caller,
                 }
            if len(good_idxs) > 0:
                t["Duplication"] = t["DTP"] / t["TP"]
            else:
                t["Duplication"] = 0
            if allow_duplicate_tp:
                sub_total = t["Total"] - t["DTP"]
            else:
                sub_total = t["Total"]
            if sub_total > 0:
                t.update({"Precision": round(float(t["TP"]) / sub_total, 4),  # Note DTP are not included
                          "Sensitivity": round(float(t["TP"]) / t["Ref"], 4),
                          "Recall": round(float(t["TP"]) / t["Ref"], 4)})
                if (t["Precision"] + t["Recall"]) > 0:
                    t.update({"F1": round(2 * ((t["Precision"] * t["Recall"]) / (t["Precision"] + t["Recall"])), 4)})
                else:
                    t["F1"] = None
            else:
                t.update({"F1": None, "Precision": None, "Recall": None})
            ts.append(t)

        else:
            df = dta[dta["strata"] >= threshold]
            if len(df) > 0:
                t = {"Total": len(df),
                     "Ref": n_in_ref,
                     "DTP": np.sum(np.in1d(df["DTP"], True)),
                     "TP": np.sum(np.in1d(df["TP"], True)),
                     "FP": np.sum(np.in1d(df["FP"], True)),
                     "T >=": threshold,
                     "Caller": data.caller,
                     }
                t["FN"] = n_in_ref - t["TP"]
                if len(good_idxs) > 0:
                    t["Duplication"] = t["DTP"] / len(good_idxs)
                else:
                    t["Duplication"] = 0
                if allow_duplicate_tp:
                    sub_total = t["Total"] - t["DTP"]
                else:
                    sub_total = t["Total"]
                if sub_total > 0:
                    t.update({"Precision": round(float(t["TP"]) / sub_total, 4),
                              "Recall": round(float(t["TP"]) / t["Ref"], 4)})
                    t.update({"F1": round(2 * ((t["Precision"] * t["Recall"]) / (t["Precision"] + t["Recall"] + 1e-6)), 4)})
                ts.append(t)

    if len(ts) == 0:
        print("Warning: precision/recall could not be determined", file=stderr)
        ts = [{"Total": None, "Ref": len(ref_bedpe), "DTP": None, "TP": None, "FP": None, "FN": None, "T >=": None,
               "Duplication": None, "Precision": None, "Recall": None, "F1": None, "strata": None}]
        data.false_negative_indexes = ref_bedpe.index

    data.scores = pd.DataFrame.from_records(ts)[["Caller", "T >=", "Ref", "Total", "TP", "FP", "DTP", "FN", "Duplication", "Precision",
                                                 "Recall", "F1"]]
    data.false_negative_indexes = missing_ref_indexes

    if show_table:

        print("Scores:", file=stderr)
        if len(data.breaks_df) == 0:
            print("None", file=stderr)
        else:
            print(data.scores.to_string(), file=stderr)

    dat = data.breaks_df
    if allow_duplicate_tp:
        dat = dat[~dat["DTP"]]
    # dat = dat[dat["quantified"]]
    if "ref_size" not in dat.columns:
        print("No TP calls found", file=stderr)
    elif "svlen" not in ref_bedpe:
        print("SV length missing from ref CallSet", file=stderr)
    else:

        s, s_ranges = np.histogram([i for i, j in zip(dat["ref_size"], dat["TP"]) if i == i and j], ref_size_bins)
        size_brackets = []
        for idx in range(len(s_ranges) - 1):
            size_brackets.append("[{}, {})".format(s_ranges[idx], s_ranges[idx + 1]))

        size_brackets.append("All ranges")
        s = np.append(s, [s.sum()])
        s2, _ = np.histogram(ref_bedpe["svlen"], ref_size_bins)
        s2 = np.append(s2, [s2.sum()])

        with np.errstate(divide='ignore', invalid='ignore'):
            sens = s / s2
            s_fp, _ = np.histogram([i for i, j in zip(dat["svlen"], dat["TP"]) if i == i and not j], ref_size_bins)
            s_fp = np.append(s_fp, [s_fp.sum()])
            prec = s / (s + s_fp)
            f1 = 2 * ((prec * sens) / (prec + sens))

        # Duplicated are currently filtered before this step
        # s_dtp, _ = np.histogram([i for i, j in zip(dat["svlen"], dat["DTP"]) if i == i and j], ref_size_bins)
        # s_dtp = np.append(s_dtp, [s_dtp.sum()])

        s_fn, _ = np.histogram([ref_data.breaks_df["svlen"].loc[i] for i in missing_ref_indexes], ref_size_bins)
        s_fn = np.append(s_fn, [s_fn.sum()])

        df_sizes = pd.DataFrame({"Caller": [data.caller] * len(size_brackets),
                                 "Ref size ranges": size_brackets, "TP": s, "FP": s_fp, "FN": s_fn, "Precision": prec,
                                 "Recall": sens, "F1": f1})
        data.size_scores = df_sizes

        if show_table:
            print("Scores over size ranges:", file=stderr)
            print(df_sizes.to_string(), file=stderr)

    if 'GT' in ref_data.breaks_df and 'GT' in data.breaks_df:
        ref_gts = ref_data.breaks_df['GT']
        gt_correct = []
        actual_gt = []
        for i, q_gt in zip(data.breaks_df.ref_index, data.breaks_df.GT):
            if i == i and i is not None:
                r_gt = ref_gts.loc[int(i)]
                if q_gt == r_gt:
                    gt_correct.append(True)
                else:
                    gt_correct.append(False)
                actual_gt.append(r_gt)
            else:
                gt_correct.append(False)
                actual_gt.append("0/0")
        data.breaks_df["ref_GT"] = actual_gt
        tp, fp = gt_correct.count(True), gt_correct.count(False)
        if tp + fp == 0:
            print("Genotype true-positives and false-positives == 0", file=stderr)

        else:
            t = {"Caller": data.caller, "GT_precision": 0 if tp + fp == 0 else tp / (tp + fp),
                 "GT_recall": 0 if len(missing_ref_indexes) == 0 and tp == 0 else tp / (tp + len(missing_ref_indexes))}

            t["GT_f1"] = round(2 * ((t["GT_precision"] * t["GT_recall"]) / (t["GT_precision"] + t["GT_recall"] + 1e-6)), 4)
            t.update({k: v for k, v in Counter(data.breaks_df.GT).items() if v == v})

            data.gt_scores = pd.DataFrame({k: [v] for k, v in t.items()})

            if show_table:
                print("GT scores:", file=stderr)
                print(data.gt_scores, file=stderr)

    # print(data.breaks_df[data.breaks_df["FP"]], file=stderr)
