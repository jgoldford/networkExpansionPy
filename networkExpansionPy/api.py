# Skeleton for a basic API using fastAPI

from fastapi import FastAPI, Query
from lib import GlobalMetabolicNetwork

app = FastAPI()


@app.get("/")
def read_root():
    return {"message": "networkExpansionPy API"}


@app.get("/expand")
def expand_network(seedSet: list[str] = Query(None)):
    """
    Expand the network given a seed set of compounds.

    Takes a list of compound IDs and returns a dictionary with two keys:
    - rxns: a list of reaction IDs
    - compounds: a list of compound IDs

    Example:
    """
    network = GlobalMetabolicNetwork()
    network.convertToIrreversible()
    network.pruneInconsistentReactions()
    # network.pruneThermodynamicallyInfeasibleReactions()
    compounds, rxns = network.expand(seedSet)
    return {"seeds": seedSet, "rxns": rxns, "compounds": compounds}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8080)
