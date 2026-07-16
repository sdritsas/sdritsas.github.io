/**
 * Casu 2.0
 * Copyright (c) 2026, Stylianos Dritsas, All Rights Reserved.
 */
export class VerticalSplitter
{
    constructor( splitterId, leftPaneId, rightPaneId, onResize )
    {
        this.splitter = document.getElementById( splitterId );
        this.leftPane = document.getElementById( leftPaneId );
        this.rightPane = document.getElementById( rightPaneId );
        this.dragging = false;
        this.startX = 0;
        this.startWidth = 0;
        this.onResize = onResize;
        this.resizeTimeout = null;
        this.setInitialRatio( );
        this.initEvents( );
    }

    setInitialRatio( )
    {
        const parent = this.splitter.parentElement;
        if ( parent && this.leftPane && this.rightPane )
        {
            const totalWidth = parent.clientWidth;
            const leftWidth = Math.round( totalWidth * 0.7 );
            this.leftPane.style.width = leftWidth + "px";
            this.rightPane.style.width = ( totalWidth - leftWidth ) + "px";
        }
    }

    initEvents( )
    {
        this.splitter.addEventListener( "mousedown", ( e ) => this.onMouseDown( e ) );
        window.addEventListener( "resize", ( ) => this.debounceResize( ) );
        this.observer = new ResizeObserver( ( ) => this.debounceResize( ) );
        this.observer.observe( this.leftPane );
    }

    onMouseDown( e )
    {
        this.dragging = true;
        this.startX = e.clientX;
        this.startWidth = this.leftPane.getBoundingClientRect( ).width;
        document.body.style.cursor = "ew-resize";
        document.addEventListener(
            "mousemove",
            ( this.onMouseMoveBound = ( ev ) => this.onMouseMove( ev ) )
        );
        document.addEventListener(
            "mouseup",
            ( this.onMouseUpBound = ( ) => this.onMouseUp( ) )
        );
    }

    onMouseMove( e )
    {
        if ( !this.dragging ) return;
        const dx = e.clientX - this.startX;
        const newWidth = Math.max( 50, this.startWidth + dx );
        this.leftPane.style.width = newWidth + "px";
        this.debounceResize( );
    }

    onMouseUp( )
    {
        this.dragging = false;
        document.body.style.cursor = "";
        document.removeEventListener( "mousemove", this.onMouseMoveBound );
        document.removeEventListener( "mouseup", this.onMouseUpBound );
        this.debounceResize( );
    }

    debounceResize( )
    {
        if ( this.resizeTimeout ) clearTimeout( this.resizeTimeout );
        this.resizeTimeout = setTimeout( ( ) =>
        {
            if ( typeof this.onResize === "function" ) this.onResize( );
        }, 50 );
    }
}