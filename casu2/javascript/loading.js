/**
 * Casu 2.0
 * Copyright (c) 2026, Stylianos Dritsas, All Rights Reserved.
 */
export class LoadingPanel {
    constructor( ) {
        this.overlay = document.getElementById( 'loading-overlay' );
        this.text = document.getElementById( 'loading-text' );
        this.spinner = document.querySelector( '.spinner' );
        this.show( );
    }

    show( message = 'Loading' ) {
        this.text.textContent = message;
        this.overlay.setAttribute( 'aria-hidden', 'false' );
    }

    update(message = undefined) {
        if (message !== undefined && this.text)
            this.text.textContent = message + ' ...';
    }

    hide( timeout = 2000 ) {
        setTimeout( ( ) => {
            this.overlay.setAttribute( 'aria-hidden', 'true' );
        }, timeout );
    }
}