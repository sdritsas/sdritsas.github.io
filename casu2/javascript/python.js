/**
 * Casu 2.0
 * Copyright (c) 2026, Stylianos Dritsas, All Rights Reserved.
 */
export async function InitializePython( )
{
    try
    {
        Casu.Loading.update( `Loading pyodide` );
        // let pyodide = await loadPyodide( { indexURL: './pyodide',
        //     stdout: msg => { } } );
        let pyodide = await loadPyodide( { indexURL: './pyodide' } );

        Casu.Loading.update( `Loading micropip` );
        await pyodide.loadPackage( 'micropip', { quiet: true } );

        const micropip = pyodide.pyimport( 'micropip' );

        const wheels = [
            'numpy-2.2.5-cp313-cp313-pyodide_2025_0_wasm32.whl',
            'pillow-11.3.0-cp313-cp313-pyodide_2025_0_wasm32.whl',
            'lxml-6.0.0-cp313-cp313-pyodide_2025_0_wasm32.whl',
            'trimesh-4.8.3-py3-none-any.whl',
            'networkx-3.4.2-py3-none-any.whl',
            'pyomanifold-0.1.0-cp313-cp313-pyodide_2025_0_wasm32.whl',
            'casu-2.0.0-py3-none-any.whl',
            ];
        for( const wheel of wheels )
        {
            const name = wheel.split( '-' )[0];
            Casu.Loading.update( `Loading ${name}` );
            await micropip.install( 'pyodide/' + wheel, { quiet: false } );
        }

        await pyodide.runPythonAsync( 'from casu import Models\nresult = Models( )' );
        Casu.Models = pyodide.globals.get( 'result' )
    }
    catch( error )
    {
        console.error( 'Error:', error );
    }
}