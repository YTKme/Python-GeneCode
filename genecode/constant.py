"""
Constant
~~~~~~~~

The `constant` module provide constant value(s).
"""

# The maximum number UID (Unique IDentifier) to retrieve
RETRIEVE_MAX = 100

# Fields available for all Sequence Databases (Nucleotide, Protein, EST, GSS)
SEQUENCE_DATABASE_FIELD_LIST = [
    {
        'search_field': ['Accession'],
        'short_field': ['ACCN'],
    },
    {
        'search_field': ['All Fields'],
        'short_field': ['ALL'],
    },
    {
        'search_field': ['Author'],
        'short_field': ['AU', 'AUTH'],
    },
    {
        'search_field': ['EC/RN Number'],
        'short_field': ['ECNO'],
    },
    {
        'search_field': ['Feature Key'], # Nucleotide, Protein, GSS
        'short_field': ['FKEY'],
    },
    {
        'search_field': ['Filter'],
        'short_field': ['FILT', 'SB'],
    },
    {
        'search_field': ['Gene Name'],
        'short_field': ['GENE'],
    },
    {
        'search_field': ['Genome Project'],
        'short_field': [],
    },
    {
        'search_field': ['Issue'],
        'short_field': ['ISS'],
    },
    {
        'search_field': ['Journal'],
        'short_field': ['JOUR'],
    },
    {
        'search_field': ['Keyword'],
        'short_field': ['KYWD'],
    },
    {
        'search_field': ['Modification Date'],
        'short_field': ['MDAT'],
    },
    {
        'search_field': ['Molecular Weight'], # Protein only
        'short_field': ['MOLWT'],
    },
    {
        'search_field': ['Organism'],
        'short_field': ['ORGN'],
    },
    {
        'search_field': ['Page Number'],
        'short_field': ['PAGE'],
    },
    {
        'search_field': ['Primary Accession'],
        'short_field': ['PACC'],
    },
    {
        'search_field': ['Primary Organism'],
        'short_field': ['PORGN'],
    },
    {
        'search_field': ['Properties'],
        'short_field': ['PROP'],
    },
    {
        'search_field': ['Protein Name'],
        'short_field': ['PROT'],
    },
    {
        'search_field': ['Publication Date'],
        'short_field': ['PDAT'],
    },
    {
        'search_field': ['SeqID String'],
        'short_field': ['SQID'],
    },
    {
        'search_field': ['Sequence Length'],
        'short_field': ['SLEN'],
    },
    {
        'search_field': ['Substance Name'],
        'short_field': ['SUBS'],
    },
    {
        'search_field': ['Text Word'],
        'short_field': ['WORD'],
    },
    {
        'search_field': ['Title'],
        'short_field': ['TI', 'TITL'],
    },
    {
        'search_field': ['Volume'],
        'short_field': ['VOL'],
    },
]

# Fields available only for EST and GSS databases
EST_GSS_DATABASE_FIELD_LIST = [
    {
        'search_field': ['Clone ID'],
    },
    {
        'search_field': ['EST Name', 'GSS Name'],
    },
    {
        'search_field': ['EST ID', 'GSS ID'],
    },
    {
        'search_field': ['Library Class'], # GSS Only
    },
    {
        'search_field': ['Library Name'], # EST Only
    },
    {
        'search_field': ['Submitter Name'],
    },
]
