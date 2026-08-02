from workflow.scripts.extract_bc_hydro_supply_validation import parse_supply_page


def test_parse_supply_page_selects_fiscal_2021_values() -> None:
    stations = [
        "G ordon M. Shrum 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
        "Revelstoke 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
        "Mica 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
        "Kootenay Canal 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
        "Peace Canyon 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
        "Seven Mile 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
        "Bridge River 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
        "Other 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
    ]
    lines = [
        "header",
        "TOTAL ELECTRICITY SALES AND SOURCES OF SUPPLY",
        "Hydroelectric generation",
        *stations,
        "9 90 9 8 80 8 7 70 7 6 60 6 5 50 5",
        "Thermal generation 1 10 1 2 20 2 3 30 3 4 40 4 5 50 5",
    ]
    while len(lines) < 28:
        lines.append("filler")
    lines.extend(
        [
            "Purchases under",
            "long-term",
            "commitments 10 1 20 2 30 3 40 4 50 5",
            "Purchases under",
            "short-term",
            "commitments 10 1 20 2 30 3 40 4 50 5",
            "Electricity trade purchases 10 1 20 2 30 3 40 4 50 5",
            "Other (10) (1) (20) (2) (30) (3) (40) (4) (50) (5)",
        ]
    )
    result = parse_supply_page("\n".join(lines))
    shrum = result.loc[result["facility"].eq("Gordon M. Shrum")].iloc[0]
    assert shrum["capacity_mw"] == 3
    assert shrum["energy_gwh"] == 30
    long_term = result.loc[result["facility"].eq("long-term commitments")].iloc[0]
    assert long_term["energy_gwh"] == 30
