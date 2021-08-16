import pathlib
import tifffile
import tqdm
import warnings
import textwrap

from quant_test_v2 import(
    validate_props,
    validate_masks,
    quantify_mask,
    format_mask_table,
    write_table,
    validate_img,
    load_marker_csv,
    quantify_channel
)

class Pipeline():
    def __init__(
        self,
        mask_paths,
        img_path=None,
        marker_csv_path=None,
        output_dir='.',
        mask_props=None,
        intensity_props=None,
        img_preprocess_func=None,
        img_preprocess_kwargs=None,
        table_prefix='',
        skip=None,
        save_RAM=True
    ):
        self.mask_paths = mask_paths
        self.img_path = img_path
        self.marker_csv_path = marker_csv_path
        self.output_dir = output_dir
        self.mask_props = mask_props
        self.intensity_props = intensity_props
        self.img_preprocess_func = img_preprocess_func
        self.img_preprocess_kwargs = img_preprocess_kwargs
        self.table_prefix = table_prefix
        self.skip = skip
        self.save_RAM = save_RAM

        assert skip in [None, 'morphology', 'intensity']
        
        validate_props(mask_props)
        self.mask_shape = validate_masks(mask_paths)

    def run(self):
        if self.skip != 'morphology':
            self.run_mask()
        if self.skip != 'intensity':
            self.run_img()


    def run_mask(self):
        if self.skip == 'intensity':
            for path in self.mask_paths:
                self._run_mask(self, path)
        else:
            self._run_mask(self, self.mask_paths[0])

    def _run_mask(self, path, flat=True):
        mask_name = pathlib.Path(path).name
        mask = tifffile.imread(path, level=0, key=0)
        print(f"Quantifying mask <{mask_name}>")
        mask_table = quantify_mask(mask, mask_props=self.mask_props)
        mask_table = format_mask_table(mask_table)
        write_table(
            mask_table, self.output_dir, 
            mask_path=path, img_path=self.img_path,
            prefix=self.table_prefix, suffix='_morphology', flat=flat
        )
        print(
            'Completed.', 
            'max id:', mask_table.index.max(), 
            'number of ids:', len(mask_table.index),
            '\n'
        )
    
    def run_img(self):
        validate_props(self.intensity_props)
        img_shape = validate_img(self.img_path)
        assert self.mask_shape == img_shape[1:3], (
            f"Mask shape ({self.mask_shape}) does not match image shape ({img_shape})"
        )
        channel_names = load_marker_csv(self.marker_csv_path)
        num_channel_names = len(channel_names)
        num_channels = img_shape[0]
        message = f'''
        Number of channel names ({num_channel_names}) does not match number of
        channels ({num_channels}) in image:
        {self.img_path}
        '''
        assert num_channel_names <= num_channels, (
            textwrap.dedent(message)
        )
        if num_channel_names != num_channels:
            warnings.warn(
                textwrap.dedent(message),
                RuntimeWarning, stacklevel=2
            )
        self.channel_names = channel_names
        if self.save_RAM:
            for path in self.mask_paths:
                self.masks = {pathlib.Path(path): tifffile.imread(path, level=0, key=0)}
                self._run_img()
        else:
            self.masks = {
                pathlib.Path(path): tifffile.imread(path, level=0, key=0)
                for path in self.mask_paths
            }
            self._run_img()

        pass
    
    def _run_img(self, flat=True):
        self.intensity_tables = {}
        mask_names = [p.name for p in self.masks.keys()]
        print(f"Quantifying channel with mask {mask_names}")
        _intensity_props = self.intensity_props[:] if self.intensity_props else []
        for i in tqdm.tqdm(range(len(self.channel_names))):
            is_first_channel = i == 0
            if is_first_channel:
                intensity_props = _intensity_props + ['centroid']
            else:
                intensity_props = _intensity_props
            
            for j, (mask_path, mask) in enumerate(self.masks.items()):
                if j == 0:
                    channel_table, processed_img = quantify_channel(
                        mask, tifffile.imread(self.img_path, level=0, key=i),
                        intensity_props=intensity_props,
                        channel_name=self.channel_names[i],
                        preprocess_func=self.img_preprocess_func,
                        preprocess_func_kwargs=self.img_preprocess_kwargs,
                        return_img=True
                    )
                    if len(self.masks) == 1: del processed_img
                else:
                    channel_table = quantify_channel(
                        mask, processed_img,
                        intensity_props=intensity_props,
                        channel_name=self.channel_names[i],
                        preprocess_func=None,
                        preprocess_func_kwargs=None,
                        return_img=False
                    )
                channel_table = format_mask_table(channel_table)
                if mask_path not in self.intensity_tables:
                    self.intensity_tables[mask_path] = channel_table
                else:
                    self.intensity_tables[mask_path].join(channel_table)
        for mask_path, table in self.intensity_tables:
            write_table(
                table, self.output_dir, 
                mask_path=mask_path, img_path=self.img_path,
                prefix=self.table_prefix, suffix='_intensity', flat=flat
            )
            print(
                'Completed.',
                mask_path.name,
                '-',
                pathlib.Path(self.img_path).name,
                '\n'
                'max id:', table.index.max(), 
                'number of ids:', len(table.index),
                '\n'
            )
