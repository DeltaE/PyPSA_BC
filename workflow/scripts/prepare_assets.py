"""
prepare_assets — stage 2 of the PyPSA-BC workflow.

Runs the per-technology asset sub-pipelines in sequence. Each is also runnable
on its own (prepare_wind.main(), prepare_solar.main(), prepare_tpp.main()).
Hydro is added last, once create_hydro_assets is migrated.

    wind    fetch CWTD -> CODERS tables -> turbine dict -> wind assets
    solar   CODERS tables -> solar assets
    tpp     CODERS tables -> tpp assets   (needs base-network buses.csv)
"""

from workflow.scripts import prepare_wind, prepare_solar, prepare_tpp


def main(force_download: bool = False):
    prepare_wind.main(force_download=force_download)
    prepare_solar.main()
    prepare_tpp.main()
    # TODO (last): prepare_hydro.main()


if __name__ == "__main__":
    main()
