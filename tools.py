#! /usr/bin/env python
# -*- coding: utf-8
import pandas as pd
import numpy as np


from Bio import Entrez
from Bio import SeqIO


def histo(df, key):
    """
    Function:
        histo( df, key )
    Arguments:
        df - pandas.DataFrame
        key - the name of the column we want to histogram.
    Description:
        Counts the occurences of distinct values in a column of data.
    Returns:
        h - a python dictionary where the keys are the distince values that
            occur and h[x] is the number of occurences.
    """
    u = list(df[key])
    h = {}
    for x in u:
        if x in h:
            h[x] += 1
        else:
            h[x] = 1

    return h


def histo_all(df):
    """
    Function:
        histo_all( df )
    Arguments:
        df - pandas.DataFrame
    Description:
        Peforms a histogram of all the columns in the DataFrame.
    Returns:
        h - a dictionary where the column names are keys and the values are the
            histograms of the respective columns.
    """
    colnames = list(df.columns.values)
    h = {}
    for col in colnames:
        h[col] = histo(df, col)

    return h


def get_genage_models(fname="./genage_models.csv"):
    """
    Function:
        get_genage_models( fname )
    Arguments:
        fname - default value is "./genage_models.csv"
    Description:
        Loads in GenAge dataset.
    Returns:
        GenAge dataset as a pandas.DataFrame.
    """
    return pd.DataFrame.from_csv(fname)


def get_sequence(seqid):
    """
    Function:
        get_sequence( seqid )
    Arguments:
        seqid - a string containing an Entrez gene id.
    Description:
        Queries the Entrez utility for GenBank records using BioPython.
    Returns:
        A parse BioPython sequence record for the entrez gene id.
    """
    Entrez.email = "thomas@synpon.com"
    handle = Entrez.efetch(
        db="nucleotide",
        rettype="gb",
        retmode="text",
        id=seqid
        )

    seq_record = SeqIO.read(handle, "gb")
    handle.close()
    return parse_seq_record(seq_record)


def parse_seq_record(seq_record):
    """
    Function:
        parse_seq_record( seq_record )
    Arguments:
        seq_record - a BioPython sequence record
    Description:
        Parses the seq_record to get the data we want from the object.
    Returns:
        Python dictionary containing desired data.
    """
    return {"seq":seq_record.seq.tostring(),
        "name": seq_record.name,
        "id": seq_record.id,
        "dbxrefs":seq_record.dbxrefs,
        "description":seq_record.description}


def get_seqids(df):
    """
    Function:
        get_seqids( df )
    Arguments:
        df - pandas.DataFrame containing GenAge dataset.
    Description:
        parses sequence ids from the GenAge DataFrame.
    Returns:
        a list of sequence ids in the proper string format.
    """
    u = list(df['entrez gene id'])
    v = list(df['entrez gene id'].isnull())
    w = []
    for k in range(len(u)):
        if not v[k]:
            w.append(str(u[k]).split('.')[0])
        else:
            w.append(None)
    return w

def get_seqs(df):
    """
    Function:
        get_seqs( df )
    Arguments:
        df - dataframe containing the GenAge dataset.
    Description:
        Collect sequence ids and sequences using get_seqids and get_sequence
    Returns:
        a list of python dictionaries representing data from seq_records.
    """
    seqids = get_seqids(df)
    res = []
    for x in seqids:
        try:
            seq_record = get_sequence(x)
        except:
            seq_record = None
        res.append(seq_record)
    return res

def make_columns(df, seqs):
    """
    Function:
        make_columns( df, seqs )
    Arguments:
        df - pandas dataframe containing GenAge dataset.
        seqs - the result of querying GenBank with BioPython for entrez gene ids
    Description:
        Gets attributes from entrez to create columns of data to append to the
        GenAge dataset, most specifically sequences.
    Returns:
        df - a copy of the GenAge dataset with sequences and other information
             added where it can be found.
    """
    attribs = None
    for x in seqs:
        if x != None:
            attribs = x.keys()
            break
    if attribs == None:
        print("seqs is a list of None values")
        return

    d = {x:[] for x in attribs}
    for record in seqs:
        if record != None:
            for attrib in attribs:
                if attrib == "dbxrefs":
                    if record[attrib] == []:
                        d[attrib].append(None)
                    else:
                        d[attrib].append(",".join(record[attrib]))
                else:
                    d[attrib].append(record[attrib])
        else:
            for attrib in attribs:
                d[attrib].append(None)

    dfnew = pd.DataFrame(d, index=df.index)

    dfnew['alt name'] = dfnew['name']
    del dfnew['name']

    for x in dfnew.columns.values:
        df[x] = dfnew[x]
    return df
