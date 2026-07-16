
/**
 * Casu 2.0
 * Copyright (c) 2026, Stylianos Dritsas, All Rights Reserved.
 */
import { VerticalSplitter } from './javascript/splitter.js';
import { LoadingPanel } from './javascript/loading.js';
import { Viewport } from './javascript/viewport.js';
import { Messages } from './javascript/messages.js';
import { Attributes } from './javascript/attributes.js';
import { InitializePython } from './javascript/python.js';

window.Casu = window.Casu || {};
window.Casu.Loading = new LoadingPanel();

window.onload = () => {
    const Casu = window.Casu;
    Casu.Viewport = new Viewport('viewport', 'view');
    Casu.Attributes = new Attributes('attributes');
    Casu.Messages = new Messages('messages');

    new VerticalSplitter('vsplit', 'view', 'tool', () => Casu.Viewport.Resize());

    InitializePython( ).then( ( ) => {
        Casu.Loading.update( 'Initializing Casu' );
        Casu.Loading.hide( );
        Casu.Attributes.Initialize( );
    });

    const bind_view = ( id, code ) =>
        document.getElementById( id ).addEventListener(
            'click', event => Casu.Viewport.SetView( code ) );
    bind_view( 'view-xy',  'xy'  );
    bind_view( 'view-yz',  'yz'  );
    bind_view( 'view-zx',  'zx'  );
    bind_view( 'view-xyz', 'xyz' );

    document.body.addEventListener('dragover', event => event.preventDefault());
    document.body.addEventListener('dragleave', _ => { });
    document.body.addEventListener('drop', event => {
        event.preventDefault();
        const file = event.dataTransfer.files[0];
        if (file) {
            const reader = new FileReader();
            reader.onload = event => {
                const arrayBuffer = event.target.result;
                const byteArray = new Uint8Array(arrayBuffer);
                Casu.Attributes.UploadFile( file.name, byteArray );
            };
            reader.readAsArrayBuffer(file);
        }
    });
};
