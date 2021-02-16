using SHA

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
end

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
) where T <: Int
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
) where {T <: Int} = fetch!(current_db(), name, global_ids, ν_min, ν_max, parameters)

fetch!(
    name            :: String,
    global_ids      :: Tuple{T, T},
    ν_min           :: Number,
    ν_max           :: Number,
    parameters      :: Union{Symbol, AbstractVector{Symbol}}=:standard
) where {T <: Int} = fetch!(current_db(), name, iso_id(global_ids...),ν_min, ν_max, parameters)

fetch!(
    name            :: String,
    global_ids      :: AbstractArray{Tuple{T, T}, 1},
    ν_min           :: Number,
    ν_max           :: Number,
    parameters      :: Union{Symbol, AbstractVector{Symbol}}=:standard
) where T <: Int = fetch!(current_db(), name, iso_id([i for i in global_ids[1]], [i for i in global_ids[2]]),ν_min, ν_max, parameters)

function iso_id(db, M::T, I::T) where T <: Union{Int, AbstractVector{Int}}
    res = DBInterface.execute(db, 
        "SELECT     global_id 
        FROM        isotopologues 
        WHERE       molecule_id IN (" * join(M, ',') * ")
        AND         local_id IN (" * join(I, ',') * ")") |> DataFrame

    if (size(res)[1] == 0)
        return missing
    else 
        return res[:, :global_id]
    end
end
iso_id(M::T, I::T) where T <: Union{Int, AbstractVector{Int}} = iso_id(current_db(), M, I)

function iso_id(db::SQLite.DB, formulas::T) where T <: Union{String, AbstractVector{String}}
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
iso_id(formulas::T) where T <: Union{String, AbstractVector{String}} = iso_id(current_db(), formulas)

function isotopologue(db::SQLite.DB, global_id::Int)    
    DBInterface.execute(db, 
        "SELECT     * 
        FROM        isotopologues
        WHERE       global_id = ?", [global_id]) |> DataFrame
end
isotopologue(global_id::Int) = isotopologue(current_db(), global_id)

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

    df = CSV.File(tmp_file; header=parameters)    

    return df
end

function print_progress(total::Integer, now::Integer)
    print("Download ", now, " of ", total, " bytes\r")
end