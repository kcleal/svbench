from svbench import CallSet, Col
import numpy as np
from sys import stderr

__all__ = ["load_dysgu", "load_lumpy", "load_delly", "load_manta", "load_sniffles", "load_whamg"]


def load_dysgu(pth, dataset, caller="dysgu"):
    c = Col("FORMAT", "PROB", bins=np.arange(0, 1.01, 0.025))
    k = Col("FILTER", op="eq", thresh=None)
    print(f"Loading dysgu, stratified using {c}, keep records with {k}", file=stderr)
    return CallSet(dataset=dataset, caller=caller).load_vcf(pth, stratify=c, keep=[k])


def load_lumpy(pth, dataset, caller="lumpy"):
    c = Col("INFO", "SU", bins=range(0, 30, 1))
    print("Loading lumpy, stratified using. (All records kept)", c, file=stderr)
    return CallSet(dataset=dataset, caller=caller).load_vcf(pth, stratify=c)


def load_delly(pth, dataset, caller="delly"):
    c = Col("QUAL", bins=range(0, 1500, 50))
    k = Col("FILTER", op="eq", thresh=None)
    print(f"Loading delly, stratified using {c}, keep records with {k}", file=stderr)
    return CallSet(dataset=dataset, caller=caller).load_vcf(pth, stratify=c, keep=[k])


def load_whamg(pth, dataset, caller="whamg"):
    c = Col("INFO", "A", bins=range(0, 20, 1))
    k = Col("FILTER", op="eq", thresh=None)
    print(f"Loading whamg, stratified using {c}, keep records with {k}", file=stderr)
    return CallSet(dataset=dataset, caller=caller).load_vcf(pth, stratify=c, keep=[k])


def load_manta(pth, dataset, caller="manta"):
    c = Col("QUAL", bins=range(0, 1500, 50))
    k = Col("FILTER", op="eq", thresh=None)
    print(f"Loading manta, stratified using {c}, keep records with {k}", file=stderr)
    return CallSet(dataset=dataset, caller=caller).load_vcf(pth, stratify=c, keep=[k])


def load_sniffles(pth, dataset, caller="sniffles"):
    c = Col("INFO", "RE", bins=np.arange(0, 20, 1))
    k = Col("FILTER", op="eq", thresh=None)
    print(f"Loading sniffles, stratified using {c}, keep records with {k}", file=stderr)
    return CallSet(dataset=dataset, caller=caller).load_vcf(pth, stratify=c, keep=[k])


def get_svim(pth, dataset, caller="svim"):
    c = Col("QUAL", bins=np.arange(0, 20, 1))
    k = Col("FILTER", op="eq", thresh=None)
    print(f"Loading svim, stratified using {c}, keep records with {k}", file=stderr)
    return CallSet(dataset=dataset, caller=caller).load_vcf(pth, stratify=c, keep=[k])
