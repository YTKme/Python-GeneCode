"""
Constant
~~~~~~~~

The `constant` module provide constant value(s).
"""

# The maximum number UID (Unique IDentifier) to retrieve
RETRIEVE_MAXIMUM = 100

# Default maximum request(s) per second for Entrez Programming Utilities
DEFAULT_MAXIMUM_REQUEST = 3
DEFAULT_MAXIMUM_REQUEST_API_KEY = 10

# Entrez Unique Identifiers (UIDs) for selected databases
ENTREZ_DATABASE = [
    {
        'entrez_atabase': 'BioProject',
        'uid_common_name': 'BioProject ID',
        'eutility_database_name': 'bioproject'
    },
    {
        'entrez_atabase': 'BioSample',
        'uid_common_name': 'BioSample ID',
        'eutility_database_name': 'biosample'
    },
    {
        'entrez_atabase': 'Books',
        'uid_common_name': 'Book ID',
        'eutility_database_name': 'books'
    },
    {
        'entrez_atabase': 'Conserved Domains',
        'uid_common_name': 'PSSM-ID',
        'eutility_database_name': 'cdd'
    },
    {
        'entrez_atabase': 'dbGaP',
        'uid_common_name': 'dbGaP ID',
        'eutility_database_name': 'gap'
    },
    {
        'entrez_atabase': 'dbVar',
        'uid_common_name': 'dbVar ID',
        'eutility_database_name': 'dbvar'
    },
    {
        'entrez_atabase': 'Gene',
        'uid_common_name': 'Gene ID',
        'eutility_database_name': 'gene'
    },
    {
        'entrez_atabase': 'Genome',
        'uid_common_name': 'Genome ID',
        'eutility_database_name': 'genome'
    },
    {
        'entrez_atabase': 'GEO Datasets',
        'uid_common_name': 'GDS ID',
        'eutility_database_name': 'gds'
    },
    {
        'entrez_atabase': 'GEO Profiles',
        'uid_common_name': 'GEO ID',
        'eutility_database_name': 'geoprofiles'
    },
    {
        'entrez_atabase': 'HomoloGene',
        'uid_common_name': 'HomoloGene ID',
        'eutility_database_name': 'homologene'
    },
    {
        'entrez_atabase': 'MeSH',
        'uid_common_name': 'MeSH ID',
        'eutility_database_name': 'mesh'
    },
    {
        'entrez_atabase': 'NCBI C++ Toolkit',
        'uid_common_name': 'Toolkit ID',
        'eutility_database_name': 'toolkit'
    },
    {
        'entrez_atabase': 'NLM Catalog',
        'uid_common_name': 'NLM Catalog ID',
        'eutility_database_name': 'nlmcatalog'
    },
    {
        'entrez_atabase': 'Nucleotide',
        'uid_common_name': 'GI number',
        'eutility_database_name': 'nuccore'
    },
    {
        'entrez_atabase': 'PopSet',
        'uid_common_name': 'PopSet ID',
        'eutility_database_name': 'popset'
    },
    {
        'entrez_atabase': 'Probe',
        'uid_common_name': 'Probe ID',
        'eutility_database_name': 'probe'
    },
    {
        'entrez_atabase': 'Protein',
        'uid_common_name': 'GI number',
        'eutility_database_name': 'protein'
    },
    {
        'entrez_atabase': 'Protein Clusters',
        'uid_common_name': 'Protein Cluster ID',
        'eutility_database_name': 'proteinclusters'
    },
    {
        'entrez_atabase': 'PubChem BioAssay',
        'uid_common_name': 'AID',
        'eutility_database_name': 'pcassay'
    },
    {
        'entrez_atabase': 'PubChem Compound',
        'uid_common_name': 'CID',
        'eutility_database_name': 'pccompound'
    },
    {
        'entrez_atabase': 'PubChem Substance',
        'uid_common_name': 'SID',
        'eutility_database_name': 'pcsubstance'
    },
    {
        'entrez_atabase': 'PubMed',
        'uid_common_name': 'PMID',
        'eutility_database_name': 'pubmed'
    },
    {
        'entrez_atabase': 'PubMed Central',
        'uid_common_name': 'PMCID',
        'eutility_database_name': 'pmc'
    },
    {
        'entrez_atabase': 'SNP',
        'uid_common_name': 'rs number',
        'eutility_database_name': 'snp'
    },
    {
        'entrez_atabase': 'SRA',
        'uid_common_name': 'SRA ID',
        'eutility_database_name': 'sra'
    },
    {
        'entrez_atabase': 'Structure',
        'uid_common_name': 'MMDB-ID',
        'eutility_database_name': 'structure'
    },
    {
        'entrez_atabase': 'Taxonomy',
        'uid_common_name': 'TaxID',
        'eutility_database_name': 'taxonomy'
    },
]

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
        'search_field': ['Feature Key'],  # Nucleotide, Protein, GSS
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
        'search_field': ['Molecular Weight'],  # Protein only
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
        'search_field': ['Library Class'],  # GSS Only
    },
    {
        'search_field': ['Library Name'],  # EST Only
    },
    {
        'search_field': ['Submitter Name'],
    },
]
