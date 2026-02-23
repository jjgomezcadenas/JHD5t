#!/usr/bin/env julia

"""
ITACA Detector Tables - LaTeX Table Generator

This script reads the activity summary CSV file and generates a LaTeX document
with formatted tables for the ITACA detector radioactivity budget.

Usage:
    julia itaca_detector_tables.jl [options]

Options:
    --input FILE     Input CSV file (default: itaca_activity_summary.csv)
    --output FILE    Output LaTeX file (default: itaca_activity_table.tex)
    --help, -h       Show this help message
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DataFrames
using CSV
using Printf
using ArgParse

# ─── LaTeX Document Generation ────────────────────────────────────────────────

"""
    format_number(x; digits=4)

Format a number for LaTeX output with appropriate precision.
"""
function format_number(x; digits=4)
    if abs(x) < 0.01
        return @sprintf("%.2e", x)
    elseif abs(x) < 1.0
        return @sprintf("%.4f", x)
    elseif abs(x) < 100.0
        return @sprintf("%.2f", x)
    elseif abs(x) < 10000.0
        return @sprintf("%.1f", x)
    else
        return @sprintf("%.0f", x)
    end
end

"""
    generate_latex_table(df::DataFrame)

Generate LaTeX table content from the activity DataFrame.
"""
function generate_latex_table(df::DataFrame)
    lines = String[]

    # Table header
    push!(lines, "\\begin{table}[htbp]")
    push!(lines, "\\centering")
    push!(lines, "\\caption{ITACA detector radioactivity budget. Activities are computed from")
    push!(lines, "         material masses and specific activities for \\textsuperscript{214}Bi and")
    push!(lines, "         \\textsuperscript{208}Tl. The ICS (Inner Copper Shield) is shown for two")
    push!(lines, "         configurations: 3~cm inner layer only, and full 15~cm with self-shielding")
    push!(lines, "         correction for the outer 12~cm.}")
    push!(lines, "\\label{tab:itaca_activity}")
    push!(lines, "\\begin{tabular}{llrrr}")
    push!(lines, "\\toprule")
    push!(lines, "Component & Material & Mass (kg) & \\textsuperscript{214}Bi (mBq) & \\textsuperscript{208}Tl (mBq) \\\\")
    push!(lines, "\\midrule")

    # Data rows
    for row in eachrow(df)
        component = row.Component
        material = row.Material

        # Skip ICS_15cm in the main table body (show separately)
        if component == "ICS_15cm"
            continue
        end

        # Handle total rows specially
        if startswith(component, "TOTAL")
            continue  # We'll add these at the end
        end

        # Format component name for LaTeX
        comp_latex = replace(component, "_" => "\\_")

        # Format numbers
        mass_str = format_number(row.Mass_kg)
        bi214_str = format_number(row.Bi214_mBq)
        tl208_str = format_number(row.Tl208_mBq)

        push!(lines, "$comp_latex & $material & $mass_str & $bi214_str & $tl208_str \\\\")
    end

    # Add midrule before totals
    push!(lines, "\\midrule")

    # Add total rows
    for row in eachrow(df)
        if row.Component == "TOTAL_ICS3cm"
            mass_str = format_number(row.Mass_kg)
            bi214_str = format_number(row.Bi214_mBq)
            tl208_str = format_number(row.Tl208_mBq)
            push!(lines, "\\textbf{Total (ICS 3cm)} & --- & $mass_str & $bi214_str & $tl208_str \\\\")
        end
    end

    push!(lines, "\\bottomrule")
    push!(lines, "\\end{tabular}")
    push!(lines, "\\end{table}")

    return join(lines, "\n")
end

"""
    generate_specific_activities_table()

Generate a LaTeX table with specific activities of all materials used in calculations.
Values are from MaterialsData.jl (Bq/kg).
"""
function generate_specific_activities_table()
    lines = String[]

    # Specific activities (Bq/kg) from MaterialsData.jl
    # Material => (Bi214, Tl208)
    materials = [
        ("Cu",       1.20e-6,   1.40e-6),
        ("PTFE",     0.10e-6,   1.28e-6),
        ("HDPE",     6.20e-6,   8.00e-6),
        ("Fe316Ti",  1.90e-3,   4.00e-4),
        ("Ti",       9.30e-4,   2.20e-4),
        ("Kapton",   80.0e-6,   110.0e-6),
    ]

    push!(lines, "\\begin{table}[htbp]")
    push!(lines, "\\centering")
    push!(lines, "\\caption{Specific activities of materials used in the ITACA detector.")
    push!(lines, "         Values are from radiopurity measurements and literature data.}")
    push!(lines, "\\label{tab:specific_activities}")
    push!(lines, "\\begin{tabular}{lrr}")
    push!(lines, "\\toprule")
    push!(lines, "Material & \\textsuperscript{214}Bi (Bq/kg) & \\textsuperscript{208}Tl (Bq/kg) \\\\")
    push!(lines, "\\midrule")

    for (mat, bi214, tl208) in materials
        bi214_str = @sprintf("%.2e", bi214)
        tl208_str = @sprintf("%.2e", tl208)
        push!(lines, "$mat & $bi214_str & $tl208_str \\\\")
    end

    push!(lines, "\\bottomrule")
    push!(lines, "\\end{tabular}")
    push!(lines, "\\end{table}")

    return join(lines, "\n")
end

"""
    generate_ics_comparison_table(df::DataFrame)

Generate a separate table comparing ICS configurations.
"""
function generate_ics_comparison_table(df::DataFrame)
    lines = String[]

    # Extract ICS rows
    ics_3cm = filter(row -> row.Component == "ICS_3cm", df)
    ics_15cm = filter(row -> row.Component == "ICS_15cm", df)
    total_3cm = filter(row -> row.Component == "TOTAL_ICS3cm", df)
    total_15cm = filter(row -> row.Component == "TOTAL_ICS15cm", df)

    push!(lines, "\\begin{table}[htbp]")
    push!(lines, "\\centering")
    push!(lines, "\\caption{Comparison of ICS configurations and their impact on total activity.")
    push!(lines, "         The full 15~cm ICS includes self-shielding correction: the outer 12~cm")
    push!(lines, "         contributes with an exponential attenuation factor based on the copper")
    push!(lines, "         attenuation length \$\\lambda \\approx 2.86\$~cm at 2.5~MeV.}")
    push!(lines, "\\label{tab:ics_comparison}")
    push!(lines, "\\begin{tabular}{lrrr}")
    push!(lines, "\\toprule")
    push!(lines, "Configuration & Mass\\textsuperscript{*} (kg) & \\textsuperscript{214}Bi (mBq) & \\textsuperscript{208}Tl (mBq) \\\\")
    push!(lines, "\\midrule")

    # ICS only
    if nrow(ics_3cm) > 0
        row = first(ics_3cm)
        push!(lines, "ICS (3 cm inner) & $(format_number(row.Mass_kg)) & $(format_number(row.Bi214_mBq)) & $(format_number(row.Tl208_mBq)) \\\\")
    end
    if nrow(ics_15cm) > 0
        row = first(ics_15cm)
        push!(lines, "ICS (15 cm with shielding) & $(format_number(row.Mass_kg)) & $(format_number(row.Bi214_mBq)) & $(format_number(row.Tl208_mBq)) \\\\")
    end

    push!(lines, "\\midrule")

    # Totals
    if nrow(total_3cm) > 0
        row = first(total_3cm)
        push!(lines, "\\textbf{Total (ICS 3 cm)} & $(format_number(row.Mass_kg)) & $(format_number(row.Bi214_mBq)) & $(format_number(row.Tl208_mBq)) \\\\")
    end
    if nrow(total_15cm) > 0
        row = first(total_15cm)
        push!(lines, "\\textbf{Total (ICS 15 cm)} & $(format_number(row.Mass_kg)) & $(format_number(row.Bi214_mBq)) & $(format_number(row.Tl208_mBq)) \\\\")
    end

    push!(lines, "\\bottomrule")
    push!(lines, "\\multicolumn{4}{l}{\\footnotesize \\textsuperscript{*}Effective mass after self-shielding correction for 15~cm configuration.} \\\\")
    push!(lines, "\\end{tabular}")
    push!(lines, "\\end{table}")

    return join(lines, "\n")
end

"""
    generate_full_latex_document(df::DataFrame)

Generate a complete LaTeX document with the activity tables.
"""
function generate_full_latex_document(df::DataFrame)
    doc = """
\\documentclass[11pt,a4paper]{article}

% Packages
\\usepackage[utf8]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage{amsmath}
\\usepackage{booktabs}
\\usepackage{siunitx}
\\usepackage[margin=2.5cm]{geometry}
\\usepackage{hyperref}

% Document info
\\title{ITACA Detector Radioactivity Budget}
\\author{Generated by itaca\\_detector\\_tables.jl}
\\date{\\today}

\\begin{document}

\\maketitle

\\section{Introduction}

This document presents the radioactivity budget for the ITACA detector, computing
the expected activities from \\textsuperscript{214}Bi (2.448~MeV gamma) and
\\textsuperscript{208}Tl (2.615~MeV gamma) for each detector component.

The detector geometry consists of:
\\begin{itemize}
    \\item \\textbf{TPC}: Cylindrical fiducial volume with \\(R = 160\\)~cm, \\(L = 150\\)~cm
    \\item \\textbf{BFD}: Barrel Field region, 5~mm PTFE shell
    \\item \\textbf{ICS}: Inner Copper Shield, 3--15~cm copper barrel and endcaps
    \\item \\textbf{FC rings}: Field Cage, 10 copper thin rings (\\(R_0 = 159.4\\)~cm mean radius, \\(A = 1.2\\)~cm\\textsuperscript{2} cross-section)
    \\item \\textbf{Light-Guides}: HDPE honeycomb structure (\\(D = 230\\)~cm, 10~cm thick)
    \\item \\textbf{FAT-GEM}: HDPE disk with holes (\\(R = 160\\)~cm, 5~mm thick, 6~mm holes at 10~mm pitch)
    \\item \\textbf{SiPM boards}: Kapton boards holding SiPM arrays (\\(D = 320\\)~cm, 0.2~mm thick)
    \\item \\textbf{Cathode}: Fe316Ti wire mesh (200~\\(\\mu\\)m wire, 5~mm pitch)
    \\item \\textbf{MARS}: Ion transport system (Ti longerons, Kapton skin/cables, HDPE ribs)
\\end{itemize}

\\section{Material Specific Activities}

Table~\\ref{tab:specific_activities} lists the specific activities of all materials
used in the ITACA detector radioactivity budget calculations. These values are based
on radiopurity measurements and literature data for \\textsuperscript{214}Bi
(from the \\textsuperscript{238}U decay chain) and \\textsuperscript{208}Tl
(from the \\textsuperscript{232}Th decay chain).

$(generate_specific_activities_table())

\\section{Activity Budget}

Table~\\ref{tab:itaca_activity} presents the activity budget for each detector component.
The total activity contribution from each component is computed as
\\(A = m \\times a\\), where \\(m\\) is the component mass and \\(a\\) is the
specific activity from Table~\\ref{tab:specific_activities}.

$(generate_latex_table(df))

\\subsection{ICS Self-Shielding}

The Inner Copper Shield (ICS) has a total thickness of 15~cm. However, gamma rays
originating from the outer layers are partially absorbed before reaching the TPC
fiducial volume. We model this self-shielding effect using exponential attenuation:

\\begin{equation}
    f_{\\text{escape}}(d) = \\exp\\left(-\\frac{d}{\\lambda}\\right)
\\end{equation}

where \\(d\\) is the depth from the inner surface and \\(\\lambda \\approx 2.86\\)~cm
is the attenuation length for copper at \\(\\sim 2.5\\)~MeV.

The inner 3~cm layer contributes fully (no self-shielding), while the outer 12~cm
is weighted by the escape probability. Table~\\ref{tab:ics_comparison} compares
the two configurations.

$(generate_ics_comparison_table(df))

\\section{Summary}

The dominant contribution to the radioactivity budget comes from the Inner Copper
Shield (ICS), which accounts for approximately 80\\% of both \\textsuperscript{214}Bi
and \\textsuperscript{208}Tl activities. Including the self-shielded contribution
from the full 15~cm ICS increases the total activity by approximately 40\\%.

\\end{document}
"""
    return doc
end

# ─── Command Line Argument Parsing ────────────────────────────────────────────

function parse_commandline()
    s = ArgParseSettings(
        description = "ITACA Detector LaTeX Table Generator",
        version = "1.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--input"
            help = "Input CSV file"
            arg_type = String
            default = "itaca_activity_summary.csv"
        "--output"
            help = "Output LaTeX file"
            arg_type = String
            default = "itaca_activity_table.tex"
    end

    return parse_args(s)
end

# ─── Main Function ────────────────────────────────────────────────────────────

function main()
    args = parse_commandline()

    input_file = args["input"]
    output_file = args["output"]

    println("ITACA Detector LaTeX Table Generator")
    println("="^50)

    # Check input file exists
    if !isfile(input_file)
        error("Input file not found: $input_file")
    end

    # Read CSV
    println("Reading: $input_file")
    df = CSV.read(input_file, DataFrame)

    println("\nData summary:")
    println("  Components: $(nrow(df))")
    println("  Columns: $(names(df))")

    # Generate LaTeX document
    println("\nGenerating LaTeX document...")
    latex_doc = generate_full_latex_document(df)

    # Write output
    open(output_file, "w") do f
        write(f, latex_doc)
    end

    println("Output written to: $output_file")

    # Print compilation instructions
    println("\nTo compile the LaTeX document:")
    println("  pdflatex $output_file")

    println("\nDone!")
end

# Run the script if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
