/**
 * Casu 2.0
 * Copyright (c) 2026, Stylianos Dritsas, All Rights Reserved.
 */
import { GUI } from './lil-gui.esm.min.js';
import noUiSlider from './nouislider.min.mjs';

export class Attributes {
    constructor(context) {
        this.ui = new GUI( { title: 'Properties',
            container: document.getElementById( context ) } );
        this.parameters = this.ui.addFolder( 'Parameters' );
        this.dictionary = { };
        this.Dispatch = {
            'fold': this.AddFold,
            'case': this.AddCase,
            'list': this.AddList,
            'dump': this.AddDump,
            'dual': this.AddDual,
            'plot': this.AddPlot,
            'bool': this.AddBool,
            'real': this.AddReal,
            'size': this.AddSize,
            'data': this.AddData,
            'down': this.AddDown,
        };
    }

    Initialize( ) {
        this.Model = 0;
        const models = Casu.Models.Available( );
        this.ui.add( this, 'ModelIndex', models.toJs( ) ).name( 'Model' );
        this.ModelChanged( );
        models.destroy( );
    }

    get ModelIndex( ) {
        return this.Model;
    }

    set ModelIndex( value ) {
        this.Model = value;
        this.ModelChanged( );
    }

    ModelChanged( ) {
        Casu.Models.Activate( this.Model );
        this.parameters.destroy( );
        this.parameters = this.ui.addFolder( 'Parameters' );
        this.attributes = JSON.parse( Casu.Models.Parameters( ) );
        Casu.Viewport.Reset( );
        this.Rebuild( this.attributes, this.parameters );
        this.Solve( );
    }

    UploadFile( name, data ) {
        const proxy = Casu.Models.Upload( name, data );
        this.Update( proxy.toJs( ) );
        proxy.destroy( );
    }

    Solve( ) {
        const result = Casu.Models.Update( this.attributes );
        this.Update( result.toJs( ) );
        result.destroy( );
    }

    Update( result ) {
        //console.log( result );
        const attributes = result['Attributes'];
        for( const attribute in attributes ) { // in !!!
            const control = this.dictionary[attribute];
            control.UpdateData( attributes[attribute] );
        }
        Casu.Viewport.Reset( );
        for( const element of result['Geometry'] ) // of !!!
            Casu.Viewport.AddMesh(
                element['Points'],
                element['Planes'],
                element['Colors'] );
        Casu.Viewport.Update( );
    }

    ParseStep( value, step ) {
        return Math.round( parseFloat( value ) / step ) * step;
    }

    ClampStep( value, min, max, step ) {
        return Math.round(Math.min(Math.max(value, min), max) / step) * step;
    }

    AddFold( settings, parameters, name ) {
        const value = settings[name];
        const group = parameters.addFolder(value.Label ?? name);
        group.title = value.Hint ?? '';
        this.dictionary[name] = group;
        this.Rebuild(value, group);
        return this;
    }

    AddCase(settings, parameters, name) {
        const value = settings[name];
        const options = value.Values;
        const control = parameters.add(value, 'Value', options)
            .name( value.Label ?? name )
            .onChange( new_value => {
                if( value.Adjusts ) {
                    for( const key of value.Adjusts.split( ' ' ) ) {
                        const adjust = settings[key];
                        const control = this.dictionary[key];
                        control.style.display = adjust.Case.includes( new_value ) ? 'flex' : 'none';
                    }
                }
                this.Solve( ); } );
        if (value.Hint) control.domElement.title = value.Hint;
        this.dictionary[name] = control;
        return this;
    }

    AddSlide(host, name, min, max, step) {
        const control = document.createElement('div');
        control.classList.add('controller', 'number', 'hasSlider');

        const title = document.createElement('div');
        title.classList.add('name');
        title.textContent = name;

        const widget = document.createElement('div');
        widget.classList.add('widget');

        const container = document.createElement('div');
        container.classList.add('slider');

        const slider = document.createElement('div');
        slider.classList.add('slider-square');

        const input = document.createElement('input');
        input.type = 'text';
        input.CasuInternal = false;

        const self = this;
        input.addEventListener('change', () => {
            if (input.CasuInternal) return;
            try {
                const values = input.value.split(' ').map(
                    value => self.ClampStep(value, min, max, step));
                slider.noUiSlider.set(values);
            } catch (except) {
                input.CasuInternal = true;
                input.value = slider.noUiSlider.get();
                input.CasuInternal = false;
            }
        });

        container.appendChild(slider);

        widget.appendChild(container);
        widget.appendChild(input);

        control.appendChild(title);
        control.appendChild(widget);

        const childrenContainer = host.domElement.querySelector('.children');
        if (childrenContainer) childrenContainer.appendChild(control);

        return [slider, input, control];
    }

    AddItem(value, items, name, index) {
        const step = value.Step ?? 1;
        const [slider, input, control] = this.AddSlide(
            items, `Item ${index}`, value.Min, value.Max, step);

        noUiSlider.create(slider, {
            animate: false,
            start: [value[`Item${index}`]],
            step: step,
            connect: 'lower',
            range: {
                'min': value.Min,
                'max': value.Max
            }
        });

        const self = this;
        slider.noUiSlider.on('update', () => {
            value[`Item${index}`] = self.ParseStep(slider.noUiSlider.get(), step);
            input.value = `${value[`Item${index}`]}`;
        });

        slider.noUiSlider.on('set', () => {
            self.Solve();
        });
        control.CasuSlider = slider;
        this.dictionary[`${name}${index}`] = control;
    }

    AddList(settings, parameters, name) {
        const value = settings[name];
        const items = parameters.addFolder(value.Label ?? name);
        if (value.Hint) items.domElement.title = value.Hint;
        for (let index = 0; index < value.Count; index++)
            this.AddItem(value, items, name, index);
        this.dictionary[name] = items;
        return this;
    }

    AddDump(settings, parameters, name) {
        const value = settings[name];
        const control = parameters.add(value, 'Value')
            .name(value.Label ?? name);
        control.UpdateData = source => control.setValue(source);
        this.dictionary[name] = control;
    }

    AddDual(settings, parameters, name) {
        const value = settings[name];
        const step = value.Step ?? 1;
        const [min, max] = [value.Min, value.Max];
        const [slider, input, control] = this.AddSlide(
            parameters, value.Label ?? name, value.Min, value.Max, step);

        noUiSlider.create(slider, {
            animate: false,
            start: value.Value,
            step: step,
            behaviour: 'drag',
            connect: true,
            range: { 'min': min, 'max': max }
        });

        const self = this;
        slider.noUiSlider.on('update', () => {
            const values = slider.noUiSlider.get().map(v => self.ClampStep(v, min, max, step));
            value.Value = values;
            input.value = values.join(' ');
        });

        slider.noUiSlider.on('set', () => {
            self.Solve();
        });
        this.dictionary[name] = control;
        return this;
    }

    AddPlot(settings, parameters, name) {
        const value = settings[name];
        const control = document.createElement('div');
        control.classList.add('controller');

        const title = document.createElement('div');
        title.classList.add('name');
        title.textContent = value.Label ?? name;

        const widget = document.createElement('div');
        widget.classList.add('widget');

        const image = document.createElement('img');
        //image.src = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAABCAYAAAC/iqxnAAAADElEQVR4nGNgGGAAAACBAAFuXGsbAAAAAElFTkSuQmCC";
        image.style.width = '100%';
        image.style.imageRendering = 'pixelated';
        image.style.objectFit = 'cover';
        control.UpdateData = source => {
            image.src = 'data:image/png;base64,' + source;
        }

        widget.appendChild(image);
        control.appendChild(title);
        control.appendChild(widget);

        const childrenContainer = parameters.domElement.querySelector('.children');
        if (childrenContainer) childrenContainer.appendChild(control);
        this.dictionary[name] = control;
    }

    AddTable(data) {
        const table = document.createElement('table');
        table.classList.add('data-table');
        for (let r = 0; r < data.length; r++) {
            const row = document.createElement('tr');
            row.classList.add(r == 0 ? 'data-table-head' : 'data-table-row')
            for (let c = 0; c < data[r].length; c++) {
                const cell = document.createElement('td');
                cell.classList.add('data-table-cell')
                cell.textContent = data[r][c];
                row.appendChild(cell);
            }
            table.appendChild(row);
        }
        return table;
    }

    AddData(settings, parameters, name) {
        const value = settings[name];
        const items = parameters.addFolder(value.Label ?? name);

        const control = document.createElement('div');
        control.classList.add('controller');
        control.CasuData = value.Value;

        const self = this;
        const table = this.AddTable(value.Value);
        control.UpdateData = source => {
            control.removeChild(control.lastChild);
            const new_table = self.AddTable(source);
            control.appendChild(new_table);
            control.CasuData = source;
        }

        control.appendChild(table);
        const childrenContainer = items.domElement.querySelector('.children');
        if (childrenContainer) childrenContainer.appendChild(control);

        items.add({
            'Download': () => {
                const csv = control.CasuData.map(row =>
                    row.map(item => `"${item}"`).join(',')).join('\n');
                const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
                const link = document.createElement('a');
                link.href = URL.createObjectURL(blob);
                link.download = 'data.csv';
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
            }
        }, 'Download');
        this.dictionary[name] = control;
    }

    AddDown( settings, parameters, name )
    {
        const value = settings[name];
        const control = parameters.add({
            'download': ( ) => {
                const result = Casu.Models.Callback( name, { } );
                const data = result.toJs( ); result.destroy( );
                const blob = new Blob([data['Data']], { type: 'application/octet-stream' });
                const link = document.createElement('a');
                link.href = URL.createObjectURL(blob);
                link.download = data['Name'];
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
            }
        }, 'download' ).name( value.Label ?? name );
        this.dictionary[name] = control;
    }

    AddBool(settings, parameters, name) {
        const value = settings[name];
        const control = parameters.add(value, 'Value')
            .name(value.Label ?? name)
            .onChange(_ => this.Solve());
        if (value.Hint) control.domElement.title = value.Hint;
        this.dictionary[name] = control;
        return this;
    }

    AddReal(settings, parameters, name) {
        const value = settings[name];
        const step = value.Step ?? 1;
        const show = value.Show ?? true;

        const [slider, input, control] = this.AddSlide(
            parameters, value.Label ?? name, value.Min, value.Max, step);

        noUiSlider.create( slider, {
            animate: false,
            start: [value.Value],
            step: step,
            connect: 'lower',
            range: {
                'min': value.Min,
                'max': value.Max
            }
        });
        control.style.display = show ? 'flex' : 'none';

        const self = this;
        slider.noUiSlider.on('update', () => {
            value.Value = self.ParseStep( slider.noUiSlider.get( ), step );

            const factor = 1 / step;
            const steps = Math.round( value.Value * factor );
            input.value = `${steps / factor}`;

        });

        slider.noUiSlider.on('set', () => {
            self.Solve();
        });

        this.dictionary[name] = control;
        return this;
    }

    AddSize(settings, parameters, name) {
        const value = settings[name];
        const step = value.Step ?? 1;
        const [slider, input, control] = this.AddSlide(
            parameters, value.Label ?? name, value.Min, value.Max, step);

        noUiSlider.create(slider, {
            animate: false,
            start: [value.Value],
            step: step,
            connect: 'lower',
            range: {
                'min': value.Min,
                'max': value.Max
            }
        });

        const self = this;
        slider.noUiSlider.on('set', () => {
            self.Solve();
        });

        slider.noUiSlider.on('update', ( ) => {
            const count = Math.round( parseFloat( slider.noUiSlider.get( ) ) / step ) * step;
            input.value = count;

            const adjust_keys = value.Adjusts.split( ' ' );
            for( const key of adjust_keys ) {
                const adjusts = settings[key];
                if( adjusts.Count !== count ) {
                    //-- Delete Existing Controls and Values
                    //--
                    for (let index = 0; index < adjusts.Count; index++) {
                        const item = `${key}${index}`;
                        if (self.dictionary[item]) {
                            self.dictionary[item].CasuSlider.noUiSlider.off();
                            delete adjusts[`Item${index}`];

                            self.dictionary[item].parentNode.removeChild(self.dictionary[item]);
                            delete self.dictionary[item];
                        }
                    }
                    //-- Create New Controls and Values
                    //--
                    const group = self.dictionary[key];
                    adjusts.Count = count;
                    settings[name].Value = count;
                    for (let index = 0; index < count; index++) {
                        const item = `${key}${index}`;
                        adjusts[`Item${index}`] = ((adjusts.Min + adjusts.Max) / 2) | 0;
                        this.AddItem(adjusts, group, key, index);
                    }
                }
            }
            // const adjusts = settings[value.Adjusts];
            // if (adjusts.Count !== count) {
            //     //-- Delete Existing Controls and Values
            //     //--
            //     for (let index = 0; index < adjusts.Count; index++) {
            //         const item = `${value.Adjusts}${index}`;
            //         if (self.dictionary[item]) {
            //             self.dictionary[item].CasuSlider.noUiSlider.off();
            //             delete adjusts[`Item${index}`];

            //             self.dictionary[item].parentNode.removeChild(self.dictionary[item]);
            //             delete self.dictionary[item];
            //         }
            //     }
            //     //-- Create New Controls and Values
            //     //--
            //     const group = self.dictionary[value.Adjusts];
            //     adjusts.Count = count;
            //     settings[name].Value = count;
            //     for (let index = 0; index < count; index++) {
            //         const item = `${value.Adjusts}${index}`;
            //         adjusts[`Item${index}`] = ((adjusts.Min + adjusts.Max) / 2) | 0;
            //         this.AddItem(adjusts, group, value.Adjusts, index);
            //     }
            // }
        });
        this.dictionary[name] = control;
        return this;
    }

    Rebuild( settings, parameters ) {
        //console.log( settings );
        for( const name in settings ) {
            const type = settings[name].Type;
            if( settings[name]['Hide'] ) continue;
            if( type in this.Dispatch ) {
                const make = this.Dispatch[type];
                make.call( this, settings, parameters, name );
            }
        }
    }
}