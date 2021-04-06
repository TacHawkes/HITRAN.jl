
mutable struct CurrentDatabase
    cur_db::Union{SQLite.DB, Nothing}
end
const CURRENT_DB = CurrentDatabase(nothing)
isdbnull() = CURRENT_DB.cur_db === nothing

function current_db()
    if isdbnull()
        @info("No custom HITRAN database specified, opening 'HITRAN.sqlite' (default)")
        open_database("HITRAN.sqlite")
    end
    return CURRENT_DB.cur_db
end
current_db(db::SQLite.DB) = (CURRENT_DB.cur_db = db)

"""
    open_database(file_path::String)

Opens the SQLite database at the given file path and sets it as the current database. 
Only use this if you want to use multiple different database.
"""
function open_database(file_path::String)
    db = SQLite.DB(file_path)

    # create default tables if this as a new database for select/filter queries           
    tables = SQLite.tables(db)
    if (haskey(tables, :name) == false || "molecules" ∉ tables[:name])        
        SQLite.load!(molecules, db, "molecules")
    end
    if (haskey(tables, :name) == false || "isotopologues" ∉ tables[:name])        
        SQLite.load!(isotopologues, db, "isotopologues")
    end
    if (haskey(tables, :name) == false || "table_hashes" ∉ tables[:name])        
        DBInterface.execute(db, 
            "CREATE TABLE table_hashes (
                table_name TEXT UNIQUE,
                query_hash TEXT
            )"
        )
    end

    current_db(db)

    return db
end

function fetch!(
    db              :: SQLite.DB,
    name            :: String,
    global_ids      :: Union{T, AbstractVector{T}},
    ν_min           :: Number,
    ν_max           :: Number,
    parameters      :: AbstractArray{String}
) where T <: Integer
    global_ids = unique(global_ids)
    
    url = build_request_url!(global_ids, ν_min, ν_max, parameters)
    # check if this is a duplicate request
    if (name in SQLite.tables(db)[:name])
        result = DBInterface.execute(db,
            "SELECT query_hash
            FROM    table_hashes
            WHERE   table_name=?",
            [name]
        ) |> DataFrame
        if size(result)[1] >= 1
            if bytes2hex(sha2_512(url)) == result[1,:query_hash]
                # identical query, skip this fetch command
                return
            end
        end
    end

    df = download_HITRAN(url, parameters)

    if (name in SQLite.tables(db)[:name])
        SQLite.drop!(db, name)
    end
    SQLite.load!(df, db, name)

    # update hash
    DBInterface.execute(db, "REPLACE INTO table_hashes VALUES (?, ?)", [name, bytes2hex(sha2_512(url))])
    
    nothing
end

"""
    fetch!([db,] name, global_ids, ν_min, ν_max, parameters)

Fetches new data from HITRANonline and stores it in the current database in the table
given by `name`. If the table with the given parameters already exists, no data download will be initiated.

# Arguments
- `db`: The database to use for storage (optional)
- `name`: The table name which can subsequently used as source table for the [`α`](@ref) function
- `global_ids`: The global isotopologue ids to consider. You can also provide a tuple (or an array of tuples) with `molecule_id`, `local_id` as identifiers
- `ν_min`: The minimum wavenumber in ``cm^{-1}`` to consider
- `ν_max`: The minimum wavenumber in ``cm^{-1}`` to consider
- `parameters`: A list of parameters to fetch. You can use parameter groups using Symbols as shortcuts, e.g. :standard for all HITRAN standard parameters. 
"""
fetch!(    
    name            :: String,
    global_ids      :: Union{T, AbstractVector{T}},
    ν_min           :: Number,
    ν_max           :: Number,
    parameters      :: AbstractArray{String}
) where {T <: Integer} = fetch!(current_db(), name, global_ids, ν_min, ν_max, parameters)

function fetch!(
    db              :: SQLite.DB,
    name            :: String,
    global_ids      :: Union{T, AbstractVector{T}},
    ν_min           :: Number,
    ν_max           :: Number,
    parameters      :: Union{Symbol, AbstractVector{Symbol}}=:standard
) where T <: Integer
    if isa(parameters, AbstractArray)
        parameters = merge_groups(:standard, parameters...)
    else
        parameters = merge_groups(:standard, parameters)
    end
    fetch!(db, name, global_ids, ν_min, ν_max, parameters)
end

fetch!(    
    name            :: String,
    global_ids      :: Union{T, AbstractVector{T}},
    ν_min           :: Number,
    ν_max           :: Number,
    parameters      :: Union{Symbol, AbstractVector{Symbol}}=:standard
) where {T <: Integer} = fetch!(current_db(), name, global_ids, ν_min, ν_max, parameters)

fetch!(
    name            :: String,
    global_ids      :: Tuple{T, T},
    ν_min           :: Number,
    ν_max           :: Number,
    parameters      :: Union{Symbol, AbstractVector{Symbol}}=:standard
) where {T <: Integer} = fetch!(current_db(), name, iso_id(global_ids...),ν_min, ν_max, parameters)

fetch!(
    name            :: String,
    global_ids      :: AbstractArray{Tuple{T, T}, 1},
    ν_min           :: Number,
    ν_max           :: Number,
    parameters      :: Union{Symbol, AbstractVector{Symbol}}=:standard
) where T <: Integer = fetch!(current_db(), name, iso_id([i[1] for i in global_ids], [i[2] for i in global_ids]),ν_min, ν_max, parameters)

"""
    iso_id([db::SQLite.DB,] M::T, I::T) where T <: Union{Integer, AbstractVector{Integer}}

Returns the global isotopologue IDs for the given molecule ids and local ids provided either as single value or as array for both parameters.

# Arguments
- `M`: molecule id or array of ids
- `I`: local isotopologue id or array of ids
"""
function iso_id(db::SQLite.DB, M, I)
    MI = zip(M, I)
    sql =   "SELECT     global_id 
            FROM        isotopologues 
            WHERE 	    (molecule_id, local_id)
            IN          (VALUES " * join(["(" * join("?"^length(t), ',') * ")" for t in MI], ',') * ")"      
    res = DBInterface.execute(db, 
        sql,
        [i for t in MI for i in t]) |> DataFrame    
    if (size(res)[1] == 0)
        return missing
    else 
        return res[:, :global_id]
    end
end
iso_id(M, I) = iso_id(current_db(), M, I)

"""
    iso_id([db::SQLite.DB,] formulas::T) where T <: Union{String, AbstractVector{String}}

Returns the global isotopologue IDs for the given molecule or isotopologue formulas.
"""
function iso_id(db::SQLite.DB, formulas)
    search_str = SQLite.esc_id(formulas)

    # get globals ids
    res = DBInterface.execute(db, 
        "SELECT     iso.global_id AS global_id
        FROM        isotopologues iso
        LEFT JOIN   molecules mol
        ON          (iso.molecule_id = mol.id)
        WHERE       (mol.formula IN (" * search_str * ") OR iso.formula IN (" * search_str * "))
        ORDER BY    iso.abundance DESC, iso.molecule_id, iso.local_id") |> DataFrame

    if (size(res)[1] == 0)
        return missing
    else 
        return res[:, :global_id]
    end
end
iso_id(formulas) = iso_id(current_db(), formulas)

function isotopologue(db::SQLite.DB, global_id)    
    DBInterface.execute(db, 
        "SELECT     * 
        FROM        isotopologues
        WHERE       global_id = ?", [global_id]) |> DataFrame
end
isotopologue(global_id) = isotopologue(current_db(), global_id)

query_local_db(db::SQLite.DB, sql::AbstractString, params=()) = DBInterface.execute(db, sql, params)
query_local_db(sql::AbstractString, params=()) = query_local_db(current_db(), sql, params)

##
const HITRAN_URL = "https://hitran.org/lbl/api?"

function build_request_url!(
    ids::Union{T, AbstractVector{T}},
    ν_min::Number,
    ν_max::Number,
    parameters::AbstractVector{String}
) where T <: Integer
    # global iso id and transition id should always be included
    if "trans_id" ∉ parameters
        pushfirst!(parameters, "trans_id")
    end
    if "global_iso_id" ∉ parameters
        pushfirst!(parameters, "global_iso_id")
    end

    id_string = join(ids, ',')
    par_string = join(parameters, ',')

    # build url
    url = @sprintf(
        "%siso_ids_list=%s&numin=%.2f&numax=%.2f&fixwidth=0&sep=[comma]&request_params=%s",
        HITRAN_URL,
        id_string,
        ν_min, ν_max,
        par_string
    )

    return url
end

function download_HITRAN(    
    url::String,
    parameters::AbstractVector{String};
    verbose=false)
        
    tmp_file = tempname()
    if verbose
        println("Download data from HITRAN at: ", url)
    end
    response = Downloads.request(
        url;
        output=tmp_file,
        progress=verbose ? print_progress : nothing, 
        throw=false        
    )
    #=if isa(response, RequestError)
        if verbose
            @error "Download failed"
            println(response)
        end
    end=#

    df = CSV.File(tmp_file; header=parameters, missingstring="#")    

    return df
end

function print_progress(total::Integer, now::Integer)
    print("Download ", now, " of ", total, " bytes\r")
end