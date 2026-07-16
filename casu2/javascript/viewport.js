/**
 * Casu 2.0
 * Copyright (c) 2026, Stylianos Dritsas, All Rights Reserved.
 */
import * as THREE from './three.module.js';
import { ArcballControls } from './arcball.js';

export class Viewport
{
    constructor( canvasId, parentId )
    {
        //-- DOM Elements
        //--
        this.canvas = document.getElementById( canvasId );
        this.parent = document.getElementById( parentId );

        //-- 3D Rendering Setup
        //--
        this.three = new THREE.WebGLRenderer(
        {
            antialias: true,
            alpha: true,
            canvas: this.canvas
        } );
        this.three.setPixelRatio( window.devicePixelRatio );
        this.three.domElement.style.background =
            'linear-gradient( 180deg, #808080 0%, #1f1f1f 100% )';
        this.scene = new THREE.Scene( );

        //-- Camera Setup
        //--
        this.zoom = 2.0;
        this.camera = new THREE.OrthographicCamera( -1, 1, 1, -1, 0.1, 10 );
        this.camera.position.set( 0, 0, 1 );
        this.camera.up.set( 0, 1, 0 );
        this.camera.lookAt( 0, 0, 0 );
        this.camera.zoom = 2.0;
        this.UpdateCamera( );

        //-- Lighting Setup
        //--
        const dirLight1 = new THREE.DirectionalLight( 0xffffff, 1 );
        dirLight1.position.set( 5, 5, 5 );
        this.scene.add( dirLight1 );

        const dirLight2 = new THREE.DirectionalLight( 0xbbbbbb, 1 );
        dirLight2.position.set( -5, -5, -5 );
        this.scene.add( dirLight2 );

        const ambientLight = new THREE.AmbientLight( 0x555555 );
        this.scene.add( ambientLight );

        //-- Controls Setup
        //--
        this.controls = new ArcballControls( this.camera, this.three.domElement, this.scene );
        this.controls.setGizmosVisible( false );
        this.controls.addEventListener( 'change', _ => this.Update( ) );

        this.Resize( );
        this.Update( );

        Casu.Viewport = this;
    }

    SetView( code ) {
        let z, y;
        switch (code) {
            case 'xy':
                z = new THREE.Vector3(0, 0, 1);
                y = new THREE.Vector3(0, 1, 0);
                break;
            case 'yz':
                z = new THREE.Vector3(1, 0, 0);
                y = new THREE.Vector3(0, 0, 1);
                break;
            case 'zx':
                z = new THREE.Vector3(0, 1, 0);
                y = new THREE.Vector3(1, 0, 0);
                break;
            case 'xyz':
                z = new THREE.Vector3(1, 1, 1).normalize();
                const right = new THREE.Vector3(0, 1, 0);
                y = new THREE.Vector3( 0.15, -0.35, 0.95 ).normalize( );
                break;
        }

        const t = this.controls._gizmos.position.clone( );
        const d = this.camera.position.clone( ).sub( t ).length( );
        const s = t.clone( ).addScaledVector( z, d );

        this.camera.position.copy( s );
        this.camera.up.copy( y );
        this.camera.lookAt( t );
        this.UpdateCamera( );

        this.controls.setCamera( this.camera );
        this.Update( );
    }

    AddMesh( points, planes, colors )
    {
        const geometry = new THREE.BufferGeometry( );
        geometry.setAttribute( 'position', new THREE.Float32BufferAttribute( points, 3 ) );
        geometry.setAttribute( 'color', new THREE.Float32BufferAttribute( colors, 4 ) );
        geometry.setIndex( planes );
        geometry.computeVertexNormals( );
        const material = new THREE.MeshStandardMaterial( {
            vertexColors: true, transparent: true, opacity:0.5, side: THREE.DoubleSide } );
        this.scene.add( new THREE.Mesh( geometry, material ) );
    }

    Reset( )
    {
        for ( let index = this.scene.children.length - 1; index >= 0; index-- )
        {
            const obj = this.scene.children[index];
            if ( !obj.isMesh ) continue;
            if ( obj.geometry ) obj.geometry.dispose( );
            if ( obj.material )
                if ( Array.isArray( obj.material ) )
                    obj.material.forEach( m => m.dispose( ) );
                else
                    obj.material.dispose( );
            this.scene.remove( obj );
        }
        this.Update( );
    }

    AspectRatio( )
    {
        return this.canvas.width / this.canvas.height;
    }

    UpdateCamera( )
    {
        const aspect = this.AspectRatio( );
        const zoom = this.zoom;
        if ( aspect >= 1 )
        {
            this.camera.left = -zoom * aspect / 2;
            this.camera.right = zoom * aspect / 2;
            this.camera.top = zoom / 2;
            this.camera.bottom = -zoom / 2;
        }
        else
        {
            this.camera.left = -zoom / 2;
            this.camera.right = zoom / 2;
            this.camera.top = zoom / aspect / 2;
            this.camera.bottom = -zoom / aspect / 2;
        }
        this.camera.updateProjectionMatrix( );
    }

    Update( )
    {
        if ( this.three && this.camera && this.scene )
        {
            this.three.render( this.scene, this.camera );
        }
    }

    Resize( )
    {
        if ( this.canvas && this.parent )
        {
            this.canvas.width = this.parent.clientWidth;
            this.canvas.height = this.parent.clientHeight;

            //console.log( `Canvas: ${this.canvas.width}x${this.canvas.height}` );
            if ( this.three && this.camera )
            {
                this.three.setSize( this.canvas.width, this.canvas.height, false );
                this.UpdateCamera( );
                this.Update( );
            }
        }
    }
}
