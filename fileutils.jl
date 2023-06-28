using CSV
using DataFrames
using Serialization
using JSON

function tryJSON(f::AbstractString)::Dict{String,Any}
    try
        return JSON.parsefile(abspath(f))
    catch e
        return Dict{String,Any}()
    end
end

function dser_path(f::AbstractString)
    return deserialize(abspath(f))
end

function ser_path(f::AbstractString,obj::Any)
    serialize(abspath(f), obj)
    return nothing
end

function read_df(f::AbstractString; kwargs...)
    return CSV.read(abspath(f), DataFrame; kwargs...)
end

function write_df(f::AbstractString, df)
    CSV.write(abspath(f),df)
end
