"""
This module provides functions to interact with the NMDC API,
including retrieving metadata into dataframes and linking across schema objects.

The :class:`NmdcClient.find` and :class:`NmdcClient.lookup` functions return data as a DataFrame.
To retrieve data as a list of dictionaries, use :class:`NmdcClient.find_dict` and
:class:`NmdcClient.lookup_dict`.
To retrieve data as a list of ``nmdc_schema.nmdc`` objects, use :class:`NmdcClient.find_full`
and :class:`NmdcClient.lookup_full`.
See the `NMDC schema documentation <https://microbiomedata.github.io/nmdc-schema/>`_
for more information on NMDC schema classes.

When working with DataFrames, there is a :class:`NmdcClient.merge_related` function
to merge related objects into one DataFrame.
This uses the more low-level :class:`NmdcClient.related_ids` function to
find related objects of a given type from a single ID.

Examples:

    These examples assume the following
    import to access Pythonic representations of the NMDC schema.

    .. code-block:: python

        >>> from nmdc_schema import nmdc

    ``nmdc_schema.nmdc`` contains the NMDC schema classes.
    Names of classes and slots can be discovered by autocompleting ``nmdc.`` in your IDE, or by
    visiting the `NMDC schema documentation <https://microbiomedata.github.io/nmdc-schema/>`_.

    To interact with the NMDC API, create an :class:`NmdcClient` instance.

    .. code-block:: python

        >>> client = NmdcClient()

    The following retrieves all studies with Wrighton as principal investigator as a DataFrame.
    Queries are specified using
    `MongoDB query syntax <https://www.mongodb.com/docs/manual/tutorial/query-documents/>`_.

    .. code-block:: python

        >>> query = {'principal_investigator.has_raw_value': {'$regex': 'Wrighton'}}
        >>> studies = client.find(nmdc.Study, query=query)
        >>> len(studies)
        2

    The following retrieves a specific biosample by its identifier and associates it
    with study metadata. The first time :class:`NmdcClient.merge_related` is called,
    it will fetch all linkages between objects in the NMDC schema and cache them in
    a file local to the package. You can force an update of this cache at a later time with
    by calling :class:`NmdcClient.fetch_links` with ``force`` set to ``True``.

    .. code-block:: python

        >>> ids = ['nmdc:bsm-13-7qxjvr77']
        >>> biosample = client.lookup(nmdc.Biosample, ids, fields=['id', 'name'])
        >>> add_study = client.merge_related(biosample, 'id', nmdc.Study, fields=['id', 'title'])
        >>> for col in add_study.columns:
        ...    print(col)
        id
        name
        id_Study
        title

    The following retrieves 10 samples then merges them with all associated WorkflowExecutions.

    .. code-block:: python

        >>> samp = client.find(nmdc.Biosample, fields=["id", "type"], limit=10)
        >>> exec = client.merge_related(samp, "id", nmdc.WorkflowExecution, fields=["id", "type"])
        >>> for col in exec.columns:
        ...     print(col)
        id
        type
        id_WorkflowExecution
        type_WorkflowExecution
"""

from enum import Enum
import json
import time
import os
import pickle
from typing import Any, cast, Dict, get_type_hints, List, Optional, Type, Tuple, TypeVar

import requests
import pandas as pd
from nmdc_schema import nmdc
from nmdc_schema.nmdc_data import get_nmdc_jsonschema_dict, get_nmdc_schema_definition
from linkml_runtime.utils.schemaview import SchemaView

SchemaClass = TypeVar("SchemaClass", bound=nmdc.NamedThing)
schema_dict = get_nmdc_jsonschema_dict()
schema_view = SchemaView(get_nmdc_schema_definition())


class SearchDirection(Enum):
    """
    Enumeration to specify the direction to follow linkages among objects
    in the schema. See :func:`NmdcClient.fetch_links` for details.
    """

    BACK = 0
    FORWARD = 1


class NmdcClient:
    """
    Class to manage calls to the NMDC API service.

    :param base_url: The base URL of the NMDC API.
    :param max_page_size: The maximum number of objects to retrieve per request.
    :param sleep_seconds: The number of seconds to wait between API requests when paginating.
    :param timeout: The number of seconds to wait for an API request to complete before
        throwing an exception.
    :param verbose: Whether to print URLs of API requests.
    :param links: A dictionary of links between objects in the NMDC schema.
        Defaults to None which fetches and caches the links on the first call
        to :func:`related_ids` or :func:`merge_related`.

    Example:

    .. code-block:: python

        >>> client = NmdcClient(max_page_size=100)
    """

    def __init__(
        self,
        base_url: str = "https://api.microbiomedata.org",
        max_page_size: int = 1000,
        sleep_seconds: float = 0.5,
        timeout: float = 30,
        verbose: bool = False,
        links: Optional[Dict[str, Tuple[List[str], List[str]]]] = None,
    ):
        self.base_url = base_url
        self.max_page_size = max_page_size
        self.sleep_seconds = sleep_seconds
        self.timeout = timeout
        self.verbose = verbose
        self.links = links
        self.codes = None

    def fetch_links(self, force: bool = False):
        """
        Construct and cache "major" links between objects in the NMDC schema.

        This function constructs links between objects in the NMDC schema
        and caches them in a pickle file.
        If the links have already been constructed and cached, they are
        loaded from the pickle file.

        The "major" links through the NMDC schema objects follow this
        directional flow. The `(process)*` notation
        represents a process that may be repeated zero or more times.

            Study (→ Study)* →
            Biosample (→ MaterialProcessing → ProcessedSample)* →
            DataGeneration → DataObject →
            (WorkflowExecution → DataObject)*

        :param force: Whether to force the construction of links even
            if they have already been cached.
        """
        if not force:
            try:
                package_dir = os.path.dirname(__file__)
                pickle_path = os.path.join(package_dir, "links.pickle")
                with open(pickle_path, "rb") as f:
                    self.links = pickle.load(f)
                    return
            except (FileNotFoundError, EOFError):
                pass

        # See https://github.com/microbiomedata/nmdc-runtime/pull/738#issuecomment-2436116400
        # for a rundown on why we look in these particular fields.

        links = {}

        def add_link(source, target):
            if source not in links:
                links[source] = ([], [])
            if target not in links:
                links[target] = ([], [])
            links[target][0].append(source)
            links[source][1].append(target)

        studies = self.find("study_set", fields=["id", "part_of"])
        for _, row in studies.iterrows():
            if isinstance(row["part_of"], list):
                for input_id in row["part_of"]:
                    add_link(input_id, row["id"])

        material_processings = self.find(
            "material_processing_set", fields=["id", "has_input", "has_output"]
        )
        for _, row in material_processings.iterrows():
            for input_id in row["has_input"]:
                add_link(input_id, row["id"])
            for output_id in row["has_output"]:
                add_link(row["id"], output_id)

        data_generations = self.find(
            "data_generation_set",
            fields=["id", "has_input", "has_output", "associated_studies"],
        )
        for _, row in data_generations.iterrows():
            for input_id in row["has_input"]:
                add_link(input_id, row["id"])
            if isinstance(row["has_output"], list):
                for output_id in row["has_output"]:
                    add_link(row["id"], output_id)
            for study_id in row["associated_studies"]:
                add_link(study_id, row["id"])

        workflow_executions = self.find(
            "workflow_execution_set",
            fields=["id", "was_informed_by", "has_input", "has_output"],
        )
        for _, row in workflow_executions.iterrows():
            for informed_id in row["was_informed_by"]:
                add_link(informed_id, row["id"])
            for input_id in row["has_input"]:
                add_link(input_id, row["id"])
            for output_id in row["has_output"]:
                add_link(row["id"], output_id)

        biosamples = self.find("biosample_set", fields=["id", "associated_studies"])
        for _, row in biosamples.iterrows():
            for study_id in row["associated_studies"]:
                add_link(study_id, row["id"])

        package_dir = os.path.dirname(__file__)
        pickle_path = os.path.join(package_dir, "links.pickle")
        with open(pickle_path, "wb") as f:
            pickle.dump(links, f)

        self.links = links

    def database_collections_for_type(self, cls: Type[SchemaClass]) -> List[str]:
        """
        Find the database collection names for a given schema class, if any.

        :param cls: The schema class to find the collection name for.
            Must be a subclass of ``nmdc_schema.nmdc.NamedThing``.
        :return: The list of collection names containing items of that class or subclass,
            or an empty list if none are found.
        """
        cls = cast(Type[nmdc.NamedThing], cls)
        descendants = [
            f"#/$defs/{name}" for name in schema_view.class_descendants(cls.class_name)
        ]
        collection = []
        for key, value in schema_dict["$defs"]["Database"]["properties"].items():
            db_items = [value["items"]]
            if "anyOf" in value["items"]:
                db_items = value["items"]["anyOf"]
            for db_item in db_items:
                if db_item["$ref"] in descendants:
                    collection.append(key)
                    break
        return collection

    def find_dict(
        self,
        cls: Type[SchemaClass] | str,
        query: Optional[Dict[str, Any]] = None,
        fields: Optional[List[str]] = None,
        limit: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        Retrieve a list of dictionaries of objects of the specified class from the NMDC schema.

        :param cls: The class type or collection of the objects to retrieve.
            Must be a subclass of ``nmdc_schema.nmdc.NamedThing`` or the name of a
            NMDC database collection.
        :param query: A dictionary representing the query to filter the results.
            Defaults to None which returns all objects.
        :param fields: A list of fields to include in the result.
            Defaults to None which returns all fields.
        :param limit: The maximum number of objects to retrieve.
            Defaults to None which retrieves all matching objects.
        :return: A list of dictionaries representing the retrieved objects.
        """

        if query is None:
            query = {}
        if fields is None:
            fields = []

        if isinstance(cls, str):
            collection = cls
        else:
            cls = cast(Type[nmdc.NamedThing], cls)
            collections = self.database_collections_for_type(cls)
            if len(collections) == 0:
                print(f"Note: No collection found for {cls}")
                return []
            if len(collections) > 1:
                print(
                    f"Note: Multiple collections found for {cls}: {collections}."
                    " Using the first one."
                )
            collection = collections[0]

            # Add the class to the filter (for heterogeneous collections)
            query["type"] = {
                "$in": [
                    schema_view.get_class(name).class_uri
                    for name in schema_view.class_descendants(cls.class_name)
                ]
            }

        max_page_size = self.max_page_size
        if limit is not None:
            max_page_size = min(limit, self.max_page_size)
        url = (
            f"{self.base_url}/nmdcschema/{collection}?"
            f"&filter={json.dumps(query)}"
            f"&projection={','.join(fields)}"
            f"&max_page_size={max_page_size}"
        )
        if self.verbose:
            print(url)
        response = requests.get(url, timeout=self.timeout)
        data = response.json()
        results = data["resources"]
        while "next_page_token" in data and (limit is None or len(results) < limit):
            page_token = data["next_page_token"]
            time.sleep(self.sleep_seconds)
            page_url = f"{url}&page_token={page_token}"
            if self.verbose:
                print(page_url)
            response = requests.get(page_url, timeout=self.timeout)
            data = response.json()
            results.extend(data["resources"])
        return results[:limit]

    def find_full(
        self,
        cls: Type[SchemaClass],
        query: Optional[Dict[str, Any]] = None,
        limit: Optional[int] = None,
    ) -> List[SchemaClass]:
        """
        Retrieve a list of full objects of the specified class from the NMDC schema.

        This function queries the NMDC API to retrieve objects of the specified class,
        converts the resulting dictionaries to instances of the specified schema class, and
        returns them as a list.

        :param cls: The class type of the objects to retrieve.
            Must be a subclass of ``nmdc_schema.nmdc.NamedThing``.
        :param query: A dictionary representing the query to filter the results.
            Defaults to None which returns all objects.
        :param limit: The maximum number of objects to retrieve.
            Defaults to None which retrieves all matching objects.
        :return: A list of instances of the specified class type.
        """
        return [
            cls(d)
            for d in self.find_dict(cls=cls, query=query, fields=None, limit=limit)
        ]

    def find(
        self,
        cls: Type[SchemaClass] | str,
        query: Optional[Dict[str, Any]] = None,
        fields: Optional[List[str]] = None,
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Retrieve a DataFrame of objects of the specified class from the NMDC schema.

        :param cls: The class type or collection of the objects to retrieve.
            Must be a subclass of ``nmdc_schema.nmdc.NamedThing`` or the name of a
            NMDC database collection.
        :param query: A dictionary representing the query to filter the objects.
            Defaults to None which returns all objects.
        :param fields: A list of fields to include in the result.
            Defaults to None which returns all fields.
        :param limit: The maximum number of objects to retrieve.
            Defaults to None which retrieves all matching objects.
        :return: A DataFrame containing the retrieved objects.
        """
        data = self.find_dict(cls=cls, query=query, fields=fields, limit=limit)
        df = pd.json_normalize(data)
        return df

    def _split_list(self, input_list: List[int], chunk_size=100):
        result = []
        for i in range(0, len(input_list), chunk_size):
            result.append(input_list[i : i + chunk_size])
        return result

    def lookup_dict(
        self,
        cls: Type[SchemaClass] | str,
        ids: List[str],
        fields: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Retrieve a list of dictionaries of objects by their identifiers.

        :param cls: The class type or collection of the objects to retrieve.
            Must be a subclass of ``nmdc_schema.nmdc.NamedThing`` or the name of a
            NMDC database collection.
        :param ids: A list of object identifiers.
        :param fields: A list of fields to include in the result.
            Defaults to None which includes all fields.
        :return: A list of dictionaries representing the retrieved objects.
        """
        chunked_ids = self._split_list(ids)
        results = []
        for chunk in chunked_ids:
            data = self.find_dict(cls=cls, query={"id": {"$in": chunk}}, fields=fields)
            results.extend(data)
        return data

    def lookup(
        self,
        cls: Type[SchemaClass] | str,
        ids: List[str],
        fields: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Retrieve a DataFrame of objects by their identifiers.

        :param cls: The class type or collection of the objects to retrieve.
            Must be a subclass of ``nmdc_schema.nmdc.NamedThing`` or the name of
            an NMDC database collection.
        :param ids: A list of object identifiers.
        :param fields: A list of fields to include in the result.
            Defaults to None which includes all fields.
        :return: A DataFrame representing the retrieved objects.
        """
        data = self.lookup_dict(cls=cls, ids=ids, fields=fields)
        df = pd.json_normalize(data)
        return df

    def lookup_full(
        self,
        cls: Type[SchemaClass] | str,
        ids: List[str],
    ) -> pd.DataFrame:
        """
        Retrieve a list of schema objects by their identifiers.

        :param cls: The class type or collection of the objects to retrieve.
            Must be a subclass of ``nmdc_schema.nmdc.NamedThing`` or the name of
            an NMDC database collection.
        :param ids: A list of object identifiers.
        :return: A list of instances of the specified class type.
        """
        data = self.lookup_dict(cls=cls, ids=ids, fields=None)
        return [cls(d) for d in data]

    def follow_links(
        self, root_id: str, direction: SearchDirection, steps: Optional[int] = None
    ) -> List[str]:
        """
        Find all object identifiers linked to the specified object identifier.

        This function finds all object identifiers linked to the specified object identifier
        by traversing the links in the specified direction.

        :param root_id: The object identifier to start from.
        :param direction: The direction to traverse the links.
            Use :class:`SearchDirection.BACK` for incoming links and
            :class:`SearchDirection.FORWARD` for outgoing links.
        :param steps: The maximum number of steps to follow.
            Defaults to None which follows all links.
        :return: A list of object identifiers linked to the specified object identifier.
        """
        if self.links is None:
            self.fetch_links()
        if root_id not in self.links:
            return []
        result = set()
        frontier = self.links[root_id][direction]
        step = 0
        while (steps is None or step < steps) and len(frontier) > 0:
            next_frontier = []
            for item_id in frontier:
                result.add(item_id)
                next_frontier.extend(self.links[item_id][direction])
            frontier = next_frontier
            step += 1
        return list(result)

    def fetch_codes(self):
        """
        Fetch and cache typecodes from the NMDC typecode API.
        """
        url = f"{self.base_url}/nmdcschema/typecodes"
        if self.verbose:
            print(url)
        self.codes = {
            d["schema_class"].split(":")[1]: d["name"]
            for d in requests.get(url, timeout=self.timeout).json()
        }

    def related_ids(
        self, ids: List[str], cls: Type[SchemaClass]
    ) -> List[Tuple[str, str]]:
        """
        Find related object identifiers of the specified class for the specified IDs.

        This function finds related objects of the specified class for the specified IDs,
        and returns a dictionary mapping each input object ID to a list of related object IDs
        matching the specified class.

        :param ids: A list of object IDs.
        :param cls: The type of objects to find related objects for.
        :return: A list of tuples mapping each source ID to a related target ID.

        Example:

        .. code-block:: python

            >>> client = NmdcClient()
            >>> related = client.related_ids(['nmdc:bsm-13-7qxjvr77'], nmdc.WorkflowExecution)
            >>> len(related)
            6
        """
        cls = cast(Type[nmdc.NamedThing], cls)
        descendants = schema_view.class_descendants(cls.class_name)
        if self.codes is None:
            self.fetch_codes()
        descendant_codes = [
            self.codes.get(c) for c in descendants if self.codes.get(c) is not None
        ]
        if self.links is None:
            self.fetch_links()
        result = []
        for source_id in ids:
            linked = set(self.follow_links(source_id, 0))
            linked.update(self.follow_links(source_id, 1))

            # We're checking if the typecode matches any of the descendant codes
            result.extend(
                [
                    (source_id, d)
                    for d in linked
                    if any(d.startswith(f"nmdc:{code}") for code in descendant_codes)
                ]
            )
        return result

    def merge_related(
        self,
        data: pd.DataFrame,
        id_column: str,
        cls: Type[SchemaClass],
        fields: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Merge related objects of the specified class into the specified DataFrame.

        This function merges related objects of the specified class into the specified DataFrame
        by finding related objects for the specified IDs and merging them into the DataFrame.

        :param data: The DataFrame to merge related objects into.
        :param id_column: The name of the column containing the object IDs.
        :param cls: The type of objects to merge into the DataFrame.
        :param fields: A list of fields to include in the result.
            Defaults to None which includes all fields.
        :return: A DataFrame containing the merged data.
        """
        # We need the ID column to perform the merge
        if fields is not None and "id" not in fields:
            fields.append("id")

        if fields is None:
            fields = []

        related = self.related_ids(data[id_column].tolist(), cls)
        related_id_column = f"{id_column}_{cls.class_name}"
        related_df = pd.DataFrame(related, columns=[id_column, related_id_column])

        with_ids_df = pd.merge(data, related_df, how="right", on=id_column)

        related_deduplicated = list(set([d[1] for d in related]))
        related_data_df = self.lookup(cls=cls, ids=related_deduplicated, fields=fields)
        related_data_df.rename(columns={"id": related_id_column}, inplace=True)

        result = pd.merge(
            with_ids_df,
            related_data_df,
            how="left",
            on=related_id_column,
            suffixes=("", f"_{cls.class_name}"),
        )

        return result


if __name__ == "__main__":
    import doctest

    doctest.testmod()
