/**
 * Casu 2.0
 * Copyright (c) 2026, Stylianos Dritsas, All Rights Reserved.
 */
export class Messages
{
    constructor( containerId )
    {
        this.container = document.getElementById( containerId );
    }

    RemoveAll( )
    {
        this.container.replaceChildren( );
    }

    AddMessage( value )
    {
        const item = document.createElement( "div" );
        item.textContent = value;
        item.className = "message-item";
        this.container.appendChild( item );
    }

    AddImage( image )
    {

        //const image = Casu.Models.Model.RenderImage();
        //console.log( image );
        // const element = document.createElement('img');
        // element.src = 'data:image/png;base64,' + image;
        // element.style.width = '50%';
        // element.style.height = 'auto';
        // document.getElementById('messages').appendChild(element);
    }
}