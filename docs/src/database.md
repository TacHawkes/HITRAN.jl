# Creating and maintaining the local HITRAN database cache

The local HITRAN database cache is stored using the SQLite format.
You can directly use and manipulate this database in order to filter for specific properties.
In future versions some easy-to-use Julia filter functions will be added. For now you can
create a view on a table (or multiple tables) and store them as a view which can be used
by the [`α`](@ref) function as a source table.

## Basic database usage

```@docs
    open_database(file_path::String)
``` 

```@docs
    fetch!
``` 

## Parameter groups

ToDo: Explanation of possible parameter groups

## Isotopologue IDs

```@docs
    iso_id
``` 